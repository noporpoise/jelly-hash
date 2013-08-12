/*
 * jellyhash.h, by Isaac Turner, August 2013, Public Domain.
*/

#ifndef JELLY_H_
#define JELLY_H_

#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <string.h>
#include <assert.h>
#include "twang.h"

// rehash limit is 2^RBITS-1, default is 15 rehashes
#define RBITS 4
typedef uint64_t HKey;
#define HASH_NULL UINT64_MAX

typedef struct
{
  volatile uint64_t *const data;
  const uint64_t dwords; // number of uint64_t words in data
  const uint32_t l, k, keylen, hwords; // hwords is num of words of k in hash
  const uint32_t kwords, klbits; // klbits is bits in top word
  const uint32_t strt_lbits, full_lbits, last_lbits;
  const uint32_t strt_kbits, full_kbits, last_kbits;
  const uint64_t nkeys, lmask;
} JellyHash;

// Functions
// number of entries is 2^l, k is number of bits per entry
static inline void jelly_hash_alloc(JellyHash *jhash, uint32_t l, uint32_t k);
static inline void jelly_hash_dealloc(JellyHash *jhash);
static inline HKey jelly_hash_find(JellyHash *jhash, const uint64_t *key,
                                   int insert, int *inserted);
static inline void jelly_hash_get_key(JellyHash *jhash, HKey loc, uint64_t *key);

// entry is: (msb) [MAX(k-l,0):entry][RBITS:rehash][1:lock] (lsb)
// all 0 is emtpy
// if lock is 1 then entry is being written

// when len == 0, UINT64_MAX >> (64-(len)) gives ~0 instead of 0
// so need to check for length == 0
#define bitmask(len)  (!!(len) * (UINT64_MAX >> (64-(len))))

// get a word from a bit array
// word64_safe avoids buffer overflow
#define word64(d,w,o) ((d)[w] >> (o) | ((o) > 0 ? (d)[(w)+1] << (64-(o)) : 0))
#define word64_safe(d,w,o,n) \
        ((d)[w] >> (o) | ((o) > 0 && (w+1)<(n) ? (d)[(w)+1] << (64-(o)) : 0))

// A possibly faster way to combine two words with a mask
//#define maskmerge(a,b,abits) ((a & abits) | (b & ~abits))
#define maskmerge(a,b,abits) (b ^ ((a ^ b) & abits))

#define MIN2(x,y) ((x) <= (y) ? (x) : (y))

// number of entries is 2^l, k is number of bits per entry
static inline void jelly_hash_alloc(JellyHash *jhash, uint32_t l, uint32_t k)
{
  assert(l>3 && l<=64 && k>3 && l<k && RBITS<32);
  uint32_t keylen, kwords, klbits;
  uint64_t nkeys, dwords, *data;
  uint32_t strt_lbits, full_lbits, last_lbits;
  uint32_t strt_kbits, full_kbits, last_kbits;

  keylen = (k > l ? k - l : 0) + RBITS + 1;
  kwords = (k+63)/64;
  klbits = k-64*(kwords-1);
  nkeys = 1UL << l;
  dwords = (keylen * nkeys + 63) / 64;

  if((data = calloc(dwords, 8)) == NULL) {
    fprintf(stderr, "[%s:%i] Error: Out of memory\n", __FILE__, __LINE__);
    exit(-1);
  }

  // Take some l bits for addressing, spread over kmer words
  last_lbits = kwords > 1 ? 1 + (l*((k&63))-1) / k : 0;
  full_lbits = kwords > 2 ? 1 + (l*63) / k : 0;
  strt_lbits = l - full_lbits - last_lbits;

  last_kbits = kwords > 1 ? klbits-last_lbits : 0;
  full_kbits = kwords > 2 ? 64-full_lbits : 0;
  strt_kbits = (kwords > 1 ? 64 : klbits);
  strt_kbits -= MIN2(strt_lbits, strt_kbits);

  printf("alloc Jelly Hash\n");
  printf("  l: %u, k: %u, keylen: %u, kwords: %u, klbits: %u, nkeys: %u\n",
         l, k, keylen, kwords, klbits, (uint32_t)nkeys);
  printf("  lbits: (%u,%u,%u)\n", strt_lbits, full_lbits, last_lbits);
  printf("  kbits: (%u,%u,%u)\n", strt_kbits, full_kbits, last_kbits);

  JellyHash j = {.data = data, .l = l, .k = k, .keylen = keylen,
                 .kwords = kwords, .klbits = klbits,
                 .hwords = (keylen+63)/64, .dwords = dwords,
                 .strt_lbits = strt_lbits, .strt_kbits = strt_kbits,
                 .full_lbits = full_lbits, .full_kbits = full_kbits,
                 .last_lbits = last_lbits, .last_kbits = last_kbits,
                 .nkeys = nkeys};
  memcpy(jhash, &j, sizeof(JellyHash));
}

static inline void jelly_hash_dealloc(JellyHash *jhash)
{
  free((void*)jhash->data);
}

static inline void _jelly_read(JellyHash *jhash, uint64_t w, uint64_t o,
                               uint64_t *restrict entry, uint32_t lenbits)
{
  const volatile uint64_t *restrict data = jhash->data;
  const uint64_t dwords = jhash->dwords, fullwords = lenbits/64;
  uint64_t i;
  for(i = 0; i < fullwords; i++, w++) entry[i] = word64_safe(data,w,o,dwords);
  if(lenbits&63) entry[i] = word64_safe(data,w,o,dwords) & bitmask(lenbits&63);
}

static inline void _jelly_write(volatile uint64_t *restrict data,
                                uint64_t pos, uint64_t w, uint64_t o,
                                uint64_t *restrict entry, uint32_t lenbits)
{
  uint64_t lastpos, lastw, lasto;
  uint64_t midw, mido, dw;
  uint64_t oldw, neww, mask = bitmask(o);
  if(o+lenbits < 64) mask |= UINT64_MAX<<(o+lenbits);

  do {
    oldw = data[w];
    neww = maskmerge(oldw, entry[0]<<o, mask);
  }
  while(!__sync_bool_compare_and_swap(data+w, oldw, neww));

  lastpos = pos+lenbits; lastw = lastpos/64; lasto = lastpos&63;

  // Middle words, word and offset of entry
  if(o > 0) { midw = 0; mido = 64 - o; }
  else { midw = 1; mido = 0; }

  for(dw = w+1; dw < lastw; dw++, midw++)
    data[dw++] = word64(entry,midw,mido);

  // Write last word
  if(lasto > 0 && lastw > w) {
    do {
      oldw = data[lastw];
      neww = word64(entry,midw,mido);
      neww = maskmerge(neww, oldw, bitmask(lasto));
    }
    while(!__sync_bool_compare_and_swap(data+lastw, oldw, neww));
  }

  // char tmp[500];
  // uint64_t i, nwords = (lenbits+63)/64, written[nwords];
  // printf("givn: %s\n", uint64s_to_binary(entry, nwords, tmp));
  // for(i = 0; i < nwords; i++) written[i] = word64(data,w+i,o);
  // printf("wrte: %s\n", uint64s_to_binary(written, nwords, tmp));
}

static inline int _jelly_acquire_lock(volatile uint64_t *restrict data,
                                      uint64_t lockw, uint64_t locko)
{
  uint64_t oldw = data[lockw];
  if((oldw>>locko) & 1) return 0; // lock already held
  uint64_t neww = oldw | (1UL<<locko);
  return __sync_bool_compare_and_swap(data+lockw, oldw, neww);
}

static inline void _jelly_release_lock(volatile uint64_t *restrict data,
                                       uint64_t lockw, uint64_t locko)
{
  uint64_t oldw, neww;
  do {
    oldw = data[lockw];
    neww = oldw & ~(1UL<<locko);
  }
  while(!__sync_bool_compare_and_swap(data+lockw, oldw, neww));
}

static inline HKey jelly_hash_find(JellyHash *jhash, const uint64_t *key,
                                   int insert, int *inserted)
{
  const uint64_t kwords = jhash->kwords, hwords = jhash->hwords, keylen = jhash->keylen;
  const uint32_t strt_lbits = jhash->strt_lbits, strt_kbits = jhash->strt_kbits;
  const uint32_t full_lbits = jhash->full_lbits, full_kbits = jhash->full_kbits;
  const uint32_t last_lbits = jhash->last_lbits, last_kbits = jhash->last_kbits;
  uint64_t i, hash[kwords], entry[hwords], found[hwords];
  uint64_t w, o, wlen, used;

  *inserted = 0;

  // hash
  for(i = 0; i+1 < kwords; i++) hash[i] = twang_mix64(key[i]);
  hash[i] = twang_mix64_mask(key[i], bitmask(jhash->klbits));

  // get loc
  HKey loc = hash[0] & bitmask(strt_lbits); hash[0]>>=strt_lbits;
  for(i=1; i+1<kwords; i++)
    loc = (loc << full_lbits) | (hash[i] & bitmask(full_lbits)); hash[i]>>=full_lbits;
  loc = (loc << last_lbits) | (hash[i] & bitmask(last_lbits)); hash[i]>>=last_lbits;

  memset(entry, 0, hwords * 8);

  // Combine remaining k-l bits for entry
  entry[0] = 2; // Rehash = 1, lock = 0
  for(i = 0, w = 0, o = RBITS+1; i < kwords; i++) {
    wlen = (i == 0 ? strt_kbits : (i+1 == kwords ? last_kbits : full_kbits));
    used = MIN2(64-o,wlen);
    entry[w] |= hash[i]<<o;
    if(used < wlen) entry[w+1] = hash[i]>>used;
    w += (o+wlen)>63;
    o = (o+wlen)&63;
  }

  // get bit at position loc
  uint64_t lockpos, lockw, locko;
  uint64_t rehash;
  volatile uint64_t *restrict const data = jhash->data;

  // +=2 because RBITS is shifted left one bit
  for(rehash = 1; rehash < (1<<RBITS); rehash++, entry[0]+=2)
  {
    // printf("loc: %zu\n", (size_t)loc);
    lockpos = keylen * loc; lockw = lockpos/64; locko = lockpos&63;

    while(1)
    {
      // Wait for lock bit (spin lock)
      while((data[lockw]>>locko) & 1);

      _jelly_read(jhash, lockw, locko, found, keylen);

      if(found[0] & 1) continue; // locked
      else if(memcmp(found, entry, hwords*8) == 0) return loc;
      else if(found[0] & (bitmask(RBITS)<<1)) break; // not empty, not match
      else {
        // Empty slot - entry not found
        if(!insert) return HASH_NULL;

        if(_jelly_acquire_lock(data, lockw, locko))
        {
          entry[0] |= 1; // Set lock = 1

          // Got lock - check still empty
          _jelly_read(jhash, lockw, locko, found, keylen);
          if(found[0] & (bitmask(RBITS)<<1))
          {
            // Not empty - check match
            int match = memcmp(found, entry, hwords*8);
            // Release lock
            _jelly_release_lock(data, lockw, locko);
            entry[0] &= ~1UL;
            if(match == 0) return loc;
            else break; // not matching -> rehash
          }

          // Insert
          _jelly_write(data, lockpos, lockw, locko, entry, keylen);
          // Need to ensure new bases are written before we release the lock
          __sync_synchronize(data); // DEV: is this needed?
          _jelly_release_lock(data, lockw, locko);
          *inserted = 1;
          return loc;
        }
        // failed to acquire lock -> try again at same loc
      }
    }

    // rehash - set new loc
    // printf(" REHASH [%c] add: %i\n", (char)key[0], (int)((rehash+1)*(rehash+2))/2);
    loc += ((rehash+1)*(rehash+2) + 0.5f)/2.0f;
    while(loc > jhash->nkeys) loc -= jhash->nkeys;
  }

  return HASH_NULL;
}


static inline void jelly_hash_get_key(JellyHash *jhash, HKey loc, uint64_t *key)
{
  const uint32_t kwords = jhash->kwords, hwords = jhash->hwords;
  const uint32_t strt_lbits = jhash->strt_lbits, strt_kbits = jhash->strt_kbits;
  const uint32_t full_lbits = jhash->full_lbits, full_kbits = jhash->full_kbits;
  const uint32_t last_lbits = jhash->last_lbits, last_kbits = jhash->last_kbits;
  uint64_t i, w, o, kbits, lbits, found[hwords+1];
  memset(found, 0, (hwords+1)*8);

  // read
  uint64_t bitpos = jhash->keylen*loc, bitw = bitpos/64, bito = bitpos&63;
  _jelly_read(jhash, bitw, bito, found, jhash->keylen);

  // char tmp[500];
  // printf("ucde: %s\n", uint64s_to_binary(found, hwords, tmp));

  // k bits
  for(i = 0, w = 0, o = 1+RBITS; i < kwords; i++)
  {
    kbits = i == 0 ? strt_kbits : (i+1 == kwords ? last_kbits : full_kbits);
    key[i] = word64(found, w, o);
    o += kbits; w += (o > 63); o &= 63;
  }

  // Undo reprobes on loc
  uint64_t sub = 0, rehashes = (found[0]>>1) & bitmask(RBITS);
  for(i = 1; i < rehashes; i++) sub += ((i+1)*(i+2))/2;  
  while(loc < sub) loc += jhash->nkeys;
  loc -= sub;
  // printf("loc: %zu\n", (size_t)loc);

  // l bits
  uint64_t loccpy = loc;
  for(i = kwords-1; i < UINT64_MAX; i--)
  {
    lbits = (i == 0 ? strt_lbits : (i+1 == kwords ? last_lbits : full_lbits));
    key[i] = (key[i]<<lbits) | (loccpy & bitmask(lbits));
    loccpy >>= lbits;
  }

  // printf("uhsh: %s\n", uint64s_to_binary(key, kwords, tmp));

  // unhash
  for(i = 0; i+1 < kwords; i++) key[i] = twang_unmix64(key[i]);
  key[i] = twang_unmix64_mask(key[i], bitmask(jhash->klbits));
}

#endif /* JELLY_H_ */
