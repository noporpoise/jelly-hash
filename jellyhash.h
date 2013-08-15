/* jellyhash.h, by Isaac Turner, August 2013, Public Domain. */
/* <turner.isaac@gmail.com> */

#ifndef JELLY_H_
#define JELLY_H_

#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <string.h>
#include <assert.h>

typedef size_t HKey;

// rehash limit is 2^(NRBITS)-1, default is 15 rehashes
#ifndef NRBITS
  #define NRBITS 4
#endif

// define the following constants before including to change word size used
#ifndef HWord
  #define HWord size_t
  #define HWORD_MAX SIZE_MAX
  #define HWORDBITS (sizeof(size_t)*8)
#endif

#define HASH_NULL HWORD_MAX

typedef struct
{
  volatile HWord *const data;
  volatile size_t collisions[1<<NRBITS], nentries;
  const size_t dwords, nbins, binsize; // dwords is num of words in data
  const size_t l, k, keylen, hwords; // hwords is num of words for keylen
  const size_t kwords, klbits; // klbits is bits in top word
} JellyHash;

// Functions
// number of entries is 2^l * b, k is number of bits per entry
static inline void jelly_hash_alloc(JellyHash *jh, uint32_t l, uint32_t b, uint32_t k);
static inline void jelly_hash_dealloc(JellyHash *jh);
static inline HKey jelly_hash_find(JellyHash *jh, const HWord *key,
                                   int insert, int *inserted);
static inline void jelly_hash_get_key(const JellyHash *jh, HKey loc, HWord *key);
static inline void jelly_hash_print_stats(const JellyHash *jh, FILE *fh);

// entry is: (msb) [k-l:entry][NRBITS:rehash][1:lock] (lsb)
// all 0 is emtpy
// if lock is 1 then entry is being written

// when len == 0, HWORD_MAX >> (HWORDBITS-(len)) gives ~0 instead of 0
// so need to check for length == 0
#define jh_bitmask(len)  ((!!(len)) * ((HWord)(HWORD_MAX) >> ((HWORDBITS)-(len))))

// get a word from a bit array
// need to multiple by !!o incase o is zero
#define jh_getword(d,w,o) \
        (((d)[w] >> (o)) | ((!!(o)) * ((d)[(w)+1] << ((HWORDBITS)-(o)))))

#define jh_nwords(x) (((x)+(HWORDBITS)-1)/(HWORDBITS))
#define jh_fullwords(x) ((x)/(HWORDBITS))
#define jh_rembits(x) ((x)&((HWORDBITS)-1)) /* returns 0..(HWORDBITS-1) */
#define jh_lastbits(x) (jh_rembits((x)+(HWORDBITS)-1)+1) /* returns 1..HWORDBITS */

// A possibly faster way to combine two words with a mask
//#define jh_maskmerge(a,b,abits) ((a & abits) | (b & ~abits))
#define jh_maskmerge(a,b,abits) ((b) ^ (((a) ^ (b)) & (abits)))

// Hash functions
#define jh_mix(x)          fasthash_mix(x)
#define jh_unmix(x)        fasthash_unmix(x)
#define jh_mix_mask(x,m)   fasthash_mix_mask(x,m)
#define jh_unmix_mask(x,m) fasthash_unmix_mask(x,m)

// Thomas Wang's hash functions
// #include "twang.h"
// #define jh_mix(x)          twang_mix64(x)
// #define jh_unmix(x)        twang_unmix64(x)
// #define jh_mix_mask(x,m)   twang_mix64_mask(x,m)
// #define jh_unmix_mask(x,m) twang_unmix64_mask(x,m)

#if HWORD_MAX < UINT64_MAX
  #undef jh_mix
  #undef jh_unmix
  #define jh_mix(x)   jh_mix_mask(x, HWORD_MAX)
  #define jh_unmix(x) jh_unmix_mask(x, HWORD_MAX)
#endif

// Mix function from Fast-Hash
// https://code.google.com/p/fast-hash/
// I've added unmix and mix/unmix mask functions
static inline HWord fasthash_mix(HWord h)
{
  h ^= h >> 23;
  h *= 0x2127599bf4325c37ULL;
  #if HWORD_MAX > UINT32_MAX
    h ^= h >> 47;
  #endif
  return h;
}

static inline HWord fasthash_unmix(HWord h)
{
  #if HWORD_MAX > UINT32_MAX
    h ^= (h >> 47);
  #endif
  h *= 11654453480509151623ULL;
  #if HWORD_MAX > UINT32_MAX
    h ^= (h >> 23) ^ (h >> 46);
  #else
    h ^= (h >> 23);
  #endif
  return h;
}

static inline HWord fasthash_mix_mask(HWord h, HWord m)
{
  // h &= m;
  h ^= h >> 23;
  h *= 0x2127599bf4325c37ULL;
  h &= m;
  #if HWORD_MAX > UINT32_MAX
    h ^= h >> 47;
  #endif
  return h;
}

static inline HWord fasthash_unmix_mask(HWord h, HWord m)
{
  #if HWORD_MAX > UINT32_MAX
    h ^= (h >> 47);
  #endif
  h *= 11654453480509151623ULL;
  h &= m;
  #if HWORD_MAX > UINT32_MAX
    h ^= (h >> 23) ^ (h >> 46);
  #else
    h ^= (h >> 23);
  #endif
  return h;
}

// number of entries is 2^l, k is number of bits per entry
static inline void jelly_hash_alloc(JellyHash *jh, uint32_t l,
                                    uint32_t binsize, uint32_t k)
{
  assert(l>3 && l<=sizeof(HWord)*8 && k>3 && binsize>0 && l<k && NRBITS<32);
  assert(jh_bitmask(HWORDBITS) == HWORD_MAX && binsize*1UL<<l < 1UL<<k);

  size_t keylen = (k > l ? k - l : 0) + NRBITS + 1;
  size_t nbins = 1UL << l;
  size_t dwords = jh_nwords(keylen * binsize * nbins);
  HWord *data = calloc(dwords + 1, sizeof(HWord));

  if(data == NULL) {
    fprintf(stderr, "[%s:%i] Error: Out of memory\n", __FILE__, __LINE__);
    exit(-1);
  }

  JellyHash j = {.data = data, .nbins = nbins, .binsize = binsize,
                 .l = l, .k = k, .keylen = keylen,
                 .kwords = jh_nwords(k), .klbits = jh_lastbits(k),
                 .hwords = jh_nwords(keylen), .dwords = dwords,
                 .collisions = {0}, .nentries = 0};

  memcpy(jh, &j, sizeof(JellyHash));
}

static inline void jelly_hash_dealloc(JellyHash *jh)
{
  free((void*)jh->data);
}

static inline void jelly_hash_print_stats(const JellyHash *jh, FILE *fh)
{
  size_t i, j, capacity = jh->nbins * jh->binsize;
  fprintf(fh, "  l: %zu, k: %zu, keylen: %zu, nbins: %zu, binsize: %zu "
              "occupancy: %zu/%zu (%.2f%%) mem: %zu bytes\n",
         jh->l, jh->k, jh->keylen, jh->nbins, jh->binsize,
         jh->nentries, capacity, 100 * (double)jh->nentries / capacity,
         jh->dwords*sizeof(HWord));

  if(jh->nentries > 0) {
    fprintf(fh, "collisions:\n");
    for(j = (1UL<<NRBITS)-1; jh->collisions[j] == 0 && j > 0; j--);
    for(i = 0; i <= j; i++)
      fprintf(fh, "%3zu: %zu\n", i, jh->collisions[i]);
  }
}


/* For debugging */
static inline char* hword_to_binary(HWord b, char *str)
{
  size_t i;
  for(i = 0; i < HWORDBITS; i++) str[i] = '0' + !!(b & (1UL<<(HWORDBITS-1-i)));
  str[HWORDBITS] = '\0';
  return str;
}

/* For debugging */
static inline char* hwords_to_binary(const HWord *b, size_t o, size_t l, char *str)
{
  size_t i;
  for(i = 0; i < l; i++)
    hword_to_binary(jh_getword(b,l-i-1,o), str+(i*HWORDBITS));
  return str;
}


static inline void _jelly_read(const JellyHash *jh, size_t w, size_t o,
                               HWord *restrict entry, size_t lenbits)
{
  const volatile HWord *restrict data = jh->data;
  size_t i, nwords = jh_nwords(lenbits);

  for(i = 0; i < nwords; i++, w++) entry[i] = jh_getword(data,w,o);
  entry[nwords-1] &= jh_bitmask(jh_lastbits(lenbits));

  // char tmp[500];
  // printf("read: %s\n", hwords_to_binary(entry, 0, jh_nwords(lenbits), tmp));
}

static inline void _jelly_write(volatile HWord *restrict data,
                                size_t pos, size_t w, size_t o,
                                HWord *restrict entry, size_t lenbits)
{
  size_t totallen, lastw, lasto;
  size_t midw, mido, dw;
  HWord oldw, neww, firstmsk, lastmsk;

  firstmsk = jh_bitmask(o);
  if(o+lenbits < HWORDBITS) firstmsk |= HWORDBITS << (o+lenbits);

  neww = entry[0]<<o;
  do {
    oldw = data[w];
    neww = jh_maskmerge(oldw, neww, firstmsk);
  }
  while(!__sync_bool_compare_and_swap(data+w, oldw, neww));

  if(o+lenbits <= HWORDBITS) return;

  totallen = pos+lenbits;
  lastw = jh_nwords(totallen)-1;
  lasto = jh_lastbits(totallen);
  lastmsk = jh_bitmask(lasto);

  // Middle words, word and offset of entry
  midw = !o; // if offset was 0, now use word 1, otherwise use word 0
  mido = jh_rembits(HWORDBITS - o);

  for(dw = w+1; dw < lastw; dw++, midw++)
    data[dw] = jh_getword(entry,midw,mido);

  // Write last word
  neww = jh_getword(entry,midw,mido);
  do {
    oldw = data[lastw];
    neww = jh_maskmerge(neww, oldw, lastmsk);
  }
  while(!__sync_bool_compare_and_swap(data+lastw, oldw, neww));

  // char tmp[500];
  // size_t nwords = jh_nwords(lenbits);
  // printf("givn: %s\n", hwords_to_binary(entry, 0, nwords, tmp));
  // printf("wrte: %s\n", hwords_to_binary((HWord*)data+w, o, nwords, tmp));
}

static inline int _jelly_acquire_lock(volatile HWord *restrict data,
                                      size_t lockw, size_t locko)
{
  HWord oldw = data[lockw];
  if((oldw>>locko) & 1) return 0; // lock already held
  HWord neww = oldw | (1UL<<locko);
  return __sync_bool_compare_and_swap(data+lockw, oldw, neww);
}

static inline void _jelly_release_lock(volatile HWord *restrict data,
                                       size_t lockw, size_t locko)
{
  HWord oldw, neww, mask = ~(1UL<<locko);
  do {
    oldw = data[lockw];
    neww = oldw & mask;
  }
  while(!__sync_bool_compare_and_swap(data+lockw, oldw, neww));
}

static inline HKey jelly_hash_find(JellyHash *jh, const HWord *key,
                                   int insert, int *inserted)
{
  volatile HWord *restrict const data = jh->data;
  const size_t kwords = jh->kwords;
  HWord hash[kwords], entry[jh->hwords+1], found[jh->hwords+1];
  size_t i, w, o, rehashes, binend, pos, lockpos, lockw, locko;
  size_t wordbits, lrem, krem, lbits, kbits;
  HKey loc;
  const HWord emptymsk = jh_bitmask(NRBITS)<<1, lastmsk = jh_bitmask(jh->klbits);

  *inserted = 0;

  // char tmpstr[500];
  // printf("inpt: %s\n", hwords_to_binary(key, 0, kwords, tmpstr));

  // hash
  for(i = 0; i+1 < kwords; i++) hash[i] = jh_mix(key[i]);
  hash[i] = jh_mix_mask(key[i] & lastmsk, lastmsk);

  // printf("hash: %s\n", hwords_to_binary(hash, 0, kwords, tmpstr));

  for(rehashes = 0; rehashes+1 < (1UL<<NRBITS); rehashes++)
  {
    // memset(entry, 0, jh->hwords*sizeof(HWord));
    entry[0] = (rehashes+1) << 1; // collisions id = 1.., lock = 0
    loc = 0;

    // split hash into loc and remaining k-l bits for entry
    HWord word;
    w = 0, o = NRBITS+1;
    lrem = jh->l; krem = jh->k;

    for(i = 0; i < kwords; i++)
    {
      word = hash[i];
      wordbits = (i+1 == kwords ? jh->klbits : HWORDBITS);
      lbits = (lrem*wordbits+krem-1) / krem;
      kbits = wordbits-lbits;
      // printf("%i kbits:%i lbits:%i\n", (int)i, (int)kbits, (int)lbits);
      loc |= (word & jh_bitmask(lbits)) << (jh->l - lrem); // l-lrem = lused
      word >>= lbits;
      entry[w] |= word << o;
      entry[w+1] = (!!o) * (word >> (HWORDBITS-o));
      w += (o+kbits) >= HWORDBITS; o = jh_rembits(o+kbits);
      lrem -= lbits; krem -= wordbits;
    }

    // printf("etry: %s\n", hwords_to_binary(entry, 0, jh->hwords, tmpstr));

    pos = jh->binsize*loc;
    lockpos = jh->keylen * pos;

    for(binend = pos + jh->binsize; pos < binend; pos++, lockpos += jh->keylen)
    {
      lockw = jh_fullwords(lockpos);
      locko = jh_rembits(lockpos);

      while(1)
      {
        // Wait for lock bit (spin lock)
        while((data[lockw]>>locko) & 1UL);

        _jelly_read(jh, lockw, locko, found, jh->keylen);

        if(memcmp(found, entry, jh->hwords*sizeof(HWord)) == 0) return pos;
        else if((found[0] & emptymsk) * !(found[0] & 1UL)) break;// !empty,!locked,!match
        else {
          // Empty or locked slot - entry not found
          if(!insert) return HASH_NULL;

          if(_jelly_acquire_lock(data, lockw, locko))
          {
            entry[0] |= 1UL; // Set lock = 1

            // Got lock - check still empty
            _jelly_read(jh, lockw, locko, found, jh->keylen);
            if(found[0] & (jh_bitmask(NRBITS)<<1))
            {
              // Not empty - check match
              int match = memcmp(found, entry, jh->hwords*sizeof(HWord));
              // Release lock
              _jelly_release_lock(data, lockw, locko);
              if(match == 0) return pos;
              else {
                entry[0] &= ~1UL;
                break; // not matching -> collision
              }
            }

            // Insert
            _jelly_write(data, lockpos, lockw, locko, entry, jh->keylen);
            // Need to ensure new bases are written before we release the lock
            __sync_synchronize(data); // DEV: is this needed?
            _jelly_release_lock(data, lockw, locko);
            *inserted = 1;

            // Update collisions
            __sync_add_and_fetch(jh->collisions+rehashes, 1);
            __sync_add_and_fetch(&jh->nentries, 1);

            return pos;
          }
        }
      } // while(1)
    } // bucket loop

    // re-hash
    for(i = 0; i+1 < kwords; i++) hash[i] = jh_mix(hash[i]);
    hash[i] = jh_mix_mask(hash[i], lastmsk);
  }

  return HASH_NULL;
}


static inline void jelly_hash_get_key(const JellyHash *jh, HKey loc, HWord *key)
{
  const size_t kwords = jh->kwords;
  size_t bitpos, bitw, bito, i, j, w, o, nrehashes;
  size_t loccpy, wordbits, lbits, kbits, lrem, krem;
  HWord found[jh->hwords+1];

  // read
  bitpos = jh->keylen*loc;
  bitw = jh_fullwords(bitpos);
  bito = jh_rembits(bitpos);
  _jelly_read(jh, bitw, bito, found, jh->keylen);

  // char tmp[500];
  // printf("inpt: %s\n", hwords_to_binary(found, 0, kwords, tmp));

  // Get bucket number
  loccpy = loc / jh->binsize;

  // merge k and l bits
  w = 0;
  o = 1+NRBITS;
  lrem = jh->l; krem = jh->k;

  for(i = 0; i < kwords; i++)
  {
    wordbits = (i+1 == kwords ? jh->klbits : HWORDBITS);
    lbits = (lrem*wordbits+krem-1) / krem;
    kbits = wordbits-lbits;
    // printf("%i kbits:%i lbits:%i\n", (int)i, (int)kbits, (int)lbits);
    key[i] = (jh_getword(found, w, o) << lbits) | (loccpy & jh_bitmask(lbits));
    loccpy >>= lbits;
    o += kbits; w += (o >= HWORDBITS); o = jh_rembits(o);
    lrem -= lbits; krem -= wordbits;
  }

  // printf("key : %s\n", hwords_to_binary(key, 0, kwords, tmp));

  // unhash
  nrehashes = ((found[0]>>1) & jh_bitmask(NRBITS)) - 1;

  for(i = 0; i <= nrehashes; i++) {
    for(j = 0; j+1 < kwords; j++) key[j] = jh_unmix(key[j]);
    key[j] = jh_unmix_mask(key[j], jh_bitmask(jh->klbits));
  }

  // printf("uhsh: %s\n", hwords_to_binary(key, 0, kwords, tmp));
}

#endif /* JELLY_H_ */
