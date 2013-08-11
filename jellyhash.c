#include "global.h"
#include <assert.h>
#include <time.h> // time() for srand()
#include <unistd.h> // getpid() for srand()

#include "jellyhash.h"
#include "twang.h"

// TODO:
// 1) Insert [done]
// 2) rehash [done]
// 3) undo rehash [done]
// 4) test, test, test!
// 5) add support for k < l (just a bit array)
// 6) tidy up writing/insert code
// 7) rehashing?

// entry is: (msb) [MAX(k-l,0):entry][RBITS:rehash][1:lock] (lsb)
// all 0 is emtpy
// if lock is 1 then entry is being written


// WORD_MAX >> (WORD_SIZE-(length)) gives WORD_MAX instead of 0 if length is 0
// need to check for length == 0
#define bitmask(len)  ((len) ? UINT64_MAX >> (64-(len)) : (uint64_t)0)
#define word64(d,w,o) ((d)[w] >> o | (d)[(w)+1] << (64-(o)))

// A possibly faster way to combine two words with a mask
//#define maskmerge(a,b,abits) ((a & abits) | (b & ~abits))
#define maskmerge(a,b,abits) (b ^ ((a ^ b) & abits))

// number of entries is 2^l, k is number of bits per entry
void jelly_hash_alloc(JellyHash *jhash, uint32_t l, uint32_t k)
{
  assert(l>3 && l<=64 && k>3 && l<k && RBITS<32);
  uint32_t keylen, kwords;
  uint64_t nkeys, data_words, *data;
  uint64_t strt_lbits, full_lbits, last_lbits;
  uint64_t strt_kbits, full_kbits, last_kbits;

  keylen = (k > l ? k - l : 0) + RBITS + 1;
  kwords = (k+64)/64;
  nkeys = 1UL << l;
  data_words = (keylen * nkeys + 63) / 64;
  if((data = calloc(data_words, 8)) == NULL) die("Out of memory");

  // Take some l bits for addressing, spread over kmer words
  last_lbits = kwords > 1 ? 1 + (l*((k&63))-1) / k : 0;
  full_lbits = kwords > 2 ? 1 + (l*63) / k : 0;
  strt_lbits = l - full_lbits - last_lbits;

  last_kbits = kwords > 1 ? jhash->klbits-last_lbits : 0;
  full_kbits = kwords > 2 ? 64-full_lbits : 0;
  strt_kbits = (kwords > 1 ? 64 : jhash->klbits);
  strt_kbits -= MIN2(strt_lbits, strt_kbits);

  JellyHash j = {.data = data, .l = l, .k = k, .keylen = keylen,
                 .kwords = kwords, .klbits = k-64*(kwords-1),
                 .hwords = (k > l ? (k-l+63)/64 : 0),
                 .strt_lbits = strt_lbits, .strt_kbits = strt_kbits,
                 .full_lbits = full_lbits, .full_kbits = full_kbits,
                 .last_lbits = last_lbits, .last_kbits = last_kbits,
                 .nkeys = nkeys};
  memcpy(jhash, &j, sizeof(JellyHash));
}

void jelly_hash_dealloc(JellyHash *jhash)
{
  free((void*)jhash->data);
}

static inline void jread(const volatile uint64_t *restrict data, uint64_t w, uint64_t o,
                        uint64_t *restrict entry, uint32_t lenbits)
{
  uint64_t i, fullwords = lenbits/64;
  for(i = 0; i < fullwords; i++, w++) entry[i] = word64(data,w,o);
  entry[i] = word64(data,w,o) & bitmask(lenbits&63);
}

HKey jelly_hash_find(JellyHash *jhash, uint64_t *key, int insert, int *inserted)
{
  const uint64_t kwords = jhash->kwords, keylen = jhash->keylen;
  const uint64_t strt_lbits = jhash->strt_lbits, strt_kbits = jhash->strt_kbits;
  const uint64_t full_lbits = jhash->full_lbits, full_kbits = jhash->full_kbits;
  const uint64_t last_lbits = jhash->last_lbits, last_kbits = jhash->last_kbits;
  uint64_t i, hash[kwords], entry[kwords], found[kwords];
  uint64_t w, o, wlen, used;

  *inserted = 0;

  // hash
  for(i = 0; i+1 < kwords; i++) hash[i] = twang_mix64(key[i]);
  hash[i] = twang_mix64_mask(key[i], bitmask(jhash->klbits));

  // get loc
  HKey loc = key[0] & bitmask(strt_lbits); key[0]>>=strt_lbits;
  for(i=1; i+1<kwords; i++)
    loc = (loc << full_lbits) | (key[i] & bitmask(full_lbits)); key[i]>>=full_lbits;
  loc = (loc << last_lbits) | (key[i] & bitmask(last_lbits)); key[i]>>=last_lbits;

  memset(entry, 0, kwords * 8);
  memset(found, 0, kwords * 8);

  // Combine remaining k-l bits for entry
  entry[0] = 1; // Rehash = 1
  for(i = 0, w = 0, o = RBITS; i < kwords; i++) {
    wlen = (i == 0 ? strt_kbits : (i+1 == kwords ? last_kbits : full_kbits));
    used = MIN2(64-o,wlen);
    entry[w] |= hash[i]<<o;
    if(used < wlen) entry[w++] = hash[i]>>used;
    w += (o+wlen)>63;
    o = (o+wlen)&63;
  }

  // get bit at position loc
  uint64_t lockpos, lockw, locko, bitpos, bitw, bito, lastpos, lastw, lasto;
  uint64_t midw, mido;
  uint64_t rehash, oldw, neww;
  volatile uint64_t *restrict const data = jhash->data;

  for(rehash = 1; rehash < 1<<RBITS; rehash++, entry[0]++)
  {
    lockpos = keylen * loc; lockw = lockpos/64; locko = lockpos&63;
    bitpos = lockpos+1; bitw = bitpos/64; bito = bitpos&63;

    while(1)
    {
      // Wait for lock bit (spin lock)
      while(((oldw = data[lockw])>>locko)&1);

      // read pos, compare
      jread(data, bitw, bito, found, keylen-1);
      if(memcmp(found, entry, kwords) == 0) return loc;
      else if((found[0] & bitmask(RBITS)) == 0)
      {
        // Empty slot - entry not found
        if(!insert) return HASH_NULL;
        // Attempt to acquire the lock
        neww = oldw | (0x1<<locko);
        uint64_t ret = __sync_val_compare_and_swap(data+lockw, oldw, neww);
        if(ret == neww)
        {
          // Got lock - check still empty
          if((data[bitw] >> bito) & bitmask(RBITS))
          {
            // Not empty - check match
            jread(data, bitw, bito, found, keylen-1);
            if(memcmp(found, entry, kwords) == 0) return loc;
            break; // not matching -> rehash
          }
          // Insert
          // Write first word
          do {
            oldw = data[bitw];
            neww = maskmerge(oldw, entry[0]<<bito, bitmask(bito));
          }
          while(__sync_val_compare_and_swap(data+bitw, oldw, neww) != neww);

          lastpos = bitpos+keylen-1; lastw = lastpos/64; lasto = lastpos&63;

          if(bito > 0) { midw = 0; mido = 64 - bito; }
          else { midw = 1; mido = 0; }

          while((++bitw) < lastw)
            data[bitw] = word64(entry,bitw,bito);

          // Write last word
          if(lastw > bitw) {
            do {
              oldw = data[lastw];
              neww = maskmerge(entry[kwords-1],oldw,bitmask(lasto));
            }
            while(__sync_val_compare_and_swap(data+bitw, oldw, neww) != neww);
          }

          // Need to ensure new bases are written before we release the lock
          __sync_synchronize(data); // DEV: is this needed?
          // release lock
          do {
            oldw = data[lockw];
            neww = oldw & ~(0x1<<locko);
          }
          while(__sync_val_compare_and_swap(data+lockw, oldw, neww) != neww);
          *inserted = 1;
          return loc;
        }
        // try again to acquire lock
      }
      else break; // slot filled: not empty or match -> rehash
    }

    // rehash - set new loc
    // if((rehash & 3) == 0) {}

    loc += (i+1)*(i+2)/2;
    if(loc > jhash->nkeys) loc -= jhash->nkeys;
  }

  return HASH_NULL;
}


void jelly_hash_get_key(JellyHash *jhash, HKey loc, uint64_t *key)
{
  const uint64_t kwords = jhash->kwords;
  const uint64_t strt_lbits = jhash->strt_lbits, strt_kbits = jhash->strt_kbits;
  const uint64_t full_lbits = jhash->full_lbits, full_kbits = jhash->full_kbits;
  const uint64_t last_lbits = jhash->last_lbits, last_kbits = jhash->last_kbits;
  uint64_t i, w, o, kbits, lbits, entry[kwords], found[kwords+1];
  memset(found, 0, (kwords+1)*8);
  memset(entry, 0, kwords*8);

  // k bits
  uint64_t bitpos = jhash->keylen*loc+1, bitw = bitpos/64, bito = bitpos&63;
  jread(jhash->data, bitw, bito, found, jhash->keylen-1);
  // split into entry
  for(i = 0, w = 0, o = 4; i < kwords; i++)
  {
    kbits = i == 0 ? strt_kbits : (i+1 == kwords ? last_kbits : full_kbits);
    entry[i] = word64(found, w, o);
    o += kbits; w += o > 63; o &= 63;
  }

  // Undo reprobes on loc
  uint64_t sub, rehashes = found[0] & bitmask(RBITS);
  for(i = 1; i < rehashes; i++) {
    sub = (i+1)*(i+2)/2;
    while(loc < sub) loc += jhash->nkeys;
    loc -= sub;
  }

  // l bits
  for(i = 0, o = 0; i < kwords; i++)
  {
    lbits = i == 0 ? strt_lbits : (i+1 == kwords ? last_lbits : full_lbits);
    entry[i] = (entry[i]<<lbits) | (loc>>o);
    o += lbits;
  }

  // unhash
  for(i = 0; i+1 < kwords; i++) key[i] = twang_unmix64(entry[i]);
  key[i] = twang_unmix64_mask(entry[i], bitmask(jhash->klbits));
}




char* uint64_to_binary(uint64_t b, char str[65])
{
  int i;
  for(i = 0; i < 64; i++) str[i] = '0' + !!(b & (1UL<<(63-i)));
  str[64] = '\0';
  return str;
}


// 1billion twang_mix64, twang_unmix64 took 41s
// 1billion twang_mix64, twang_unmix64_alt took 46s
// 1billion twang_mix64_mask + twang_unmix64_mask took 36s

static void test_twang()
{
  uint64_t i, n = 1000000000UL, word, hash, unhash;
  uint64_t mask = 1UL<<(int)((rand()*63.0)/RAND_MAX)-1;
  for(i = 0; i < n; i++)
  {
    // word = ((uint64_t)rand()) << 32 | rand();
    word = i;
    word &= mask;
    hash = twang_mix64_mask(word, mask);
    unhash = twang_unmix64_mask(hash, mask);
    if(word != unhash) die("%zu", (size_t)i);
  }
}

int main()
{
  // Initiate rand
  srand((unsigned int)time(NULL) + getpid());
  test_twang();
  return EXIT_SUCCESS;
}
