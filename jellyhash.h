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

// rehash limit is 2^(JH_NRBITS)-1, default is 15 rehashes
#ifndef JH_NRBITS
  #define JH_NRBITS 4
#endif

// define the following constants before including to change word size used
#if defined(JH_FORCE_64)
  #define HWord uint64_t
  #define HWORD_MAX UINT64_MAX
  #define HWORDBITS 64
#elif defined(JH_FORCE_32)
  #define HWord uint32_t
  #define HWORD_MAX UINT32_MAX
  #define HWORDBITS 32
#else
  #define HWord size_t
  #define HWORD_MAX SIZE_MAX
  #define HWORDBITS (sizeof(size_t)*8)
#endif

#define JHASH_NULL HWORD_MAX

typedef struct
{
  volatile HWord *const data;
  volatile size_t collisions[1<<JH_NRBITS], nentries;
  const size_t dwords, nbins, binsize; // dwords is num of words in data
  const size_t l, k, keylen, hwords; // hwords is num of words for keylen
  const size_t kwords, klbits; // klbits is bits in top word
} JellyHash;

// Functions
// number of entries is 2^l * b, k is number of bits per entry
static inline void jelly_hash_alloc(JellyHash *jh, uint32_t l, uint32_t b, uint32_t k);
static inline void jelly_hash_dealloc(JellyHash *jh);
static inline HKey jelly_hash_find(JellyHash *jh, const char *key,
                                   int insert, int *inserted);
static inline void jelly_hash_get_key(const JellyHash *jh, HKey loc, char *key);
static inline void jelly_hash_print_stats(const JellyHash *jh, FILE *fh);

// entry is: (msb) [k-l:entry][JH_NRBITS:rehash][1:lock] (lsb)
// all 0 is emtpy
// if lock is 1 then entry is being written

// when len == 0, HWORD_MAX >> (HWORDBITS-(len)) gives ~0 instead of 0
// so need to check for length == 0
#define jh_bitmask(len)  ((!!(len)) * ((HWord)(HWORD_MAX) >> ((HWORDBITS)-(len))))

// get a word from a bit array
// need to multiple by !!o in case o is zero
#define jh_getword(d,w,o) \
        (((d)[w] >> (o)) | ((!!(o)) * ((d)[(w)+1] << ((HWORDBITS)-(o)))))

#define jh_nwords(x) (((x)+(HWORDBITS)-1)/(HWORDBITS))
#define jh_fullwords(x) ((x)/(HWORDBITS))
#define jh_rembits(x) ((x)&((HWORDBITS)-1)) /* returns 0..(HWORDBITS-1) */
#define jh_lastbits(x) (jh_rembits((x)+(HWORDBITS)-1)+1) /* returns 1..HWORDBITS */

// A possibly faster way to combine two words with a mask
//#define jh_maskmerge(a,b,abits) ((a & abits) | (b & ~abits))
#define jh_maskmerge(a,b,abits) ((b) ^ (((a) ^ (b)) & (abits)))

#define JH_HALFBITS (HWORDBITS/2)
static const HWord jh_botmsk = HWORD_MAX >> JH_HALFBITS;
static const HWord jh_topmsk = HWORD_MAX << JH_HALFBITS;

// Hash functions

#if HWORD_MAX > UINT32_MAX
  #define jh_mix(x)            fasthash64_mix(x)
  #define jh_unmix(x)          fasthash64_unmix(x)
  // #define jh_mix_mask(x,k,m)          fasthash64_mix_mask(x,m)
  // #define jh_unmix_mask(x,k,m)        fasthash64_unmix_mask(x,m)
  #define jh_mix_mask(x,k,m)   ((k) <= 32 ? twang32_mix_mask(x,m) : fasthash64_mix_mask(x,m))
  #define jh_unmix_mask(x,k,m) ((k) <= 32 ? twang32_unmix_mask(x,m) : fasthash64_unmix_mask(x,m))
#else
  #define jh_mix(x)            twang32_mix(x)
  #define jh_unmix(x)          twang32_unmix(x)
  #define jh_mix_mask(x,k,m)   twang32_mix_mask(x,m)
  #define jh_unmix_mask(x,k,m) twang32_unmix_mask(x,m)
#endif

// more reversible hash functions at:
// https://github.com/facebook/folly/blob/master/folly/Hash.h

// Mix function from Fast-Hash used for 64bit
// https://code.google.com/p/fast-hash/
// I've added unmix and mix/unmix mask functions
static inline uint64_t fasthash64_mix(uint64_t h)
{
  h ^= h >> 23;
  h *= 0x2127599bf4325c37ULL;
  h ^= h >> 47;
  return h;
}

static inline uint64_t fasthash64_unmix(uint64_t h)
{
  h ^= (h >> 47);
  h *= 11654453480509151623ULL;
  h ^= (h >> 23) ^ (h >> 46);
  return h;
}

static inline uint64_t fasthash64_mix_mask(uint64_t h, uint64_t m)
{
  h ^= h >> 23;
  h *= 0x2127599bf4325c37ULL;
  h &= m;
  h ^= h >> 47;
  return h;
}

static inline uint64_t fasthash64_unmix_mask(uint64_t h, uint64_t m)
{
  h ^= (h >> 47);
  h *= 11654453480509151623ULL;
  h &= m;
  h ^= (h >> 23) ^ (h >> 46);
  return h;
}

static inline uint32_t twang32_mix(uint32_t h)
{
  h += ~(h<<15);
  h ^=  (h>>10);
  h +=  (h<<3);
  h ^=  (h>>6);
  h += ~(h<<11);
  h ^=  (h>>16);
  return h;
}

static inline uint32_t twang32_unmix(uint32_t h)
{
  h ^= (h>>16);
  h  = (h+1) * 4196353U;
  h ^= (h>>6) ^ (h>>12) ^ (h>>18) ^ (h>>24) ^ (h>>30);
  h += h * 954437176U;
  h ^= (h>>10) ^ (h>>20) ^ (h>>30);
  h  = (h+1) * 1073774593;
  return h;
}

static inline uint32_t twang32_mix_mask(uint32_t h, uint32_t m)
{
  h += ~(h<<15);
  h &=  m;
  h ^=  (h>>10);
  h +=  (h<<3);
  h &=  m;
  h ^=  (h>>6);
  h += ~(h<<11);
  h &=  m;
  h ^=  (h>>16);
  return h;
}

static inline uint32_t twang32_unmix_mask(uint32_t h, uint32_t m)
{
  h ^= (h>>16);
  h  = (h+1) * 4196353U;
  h &= m;
  h ^= (h>>6) ^ (h>>12) ^ (h>>18) ^ (h>>24) ^ (h>>30);
  h += h * 954437176U;
  h &= m;
  h ^= (h>>10) ^ (h>>20) ^ (h>>30);
  h  = (h+1) * 1073774593;
  h &= m;
  return h;
}

// number of entries is 2^l, k is number of bits per entry
static inline void jelly_hash_alloc(JellyHash *jh, uint32_t l,
                                    uint32_t binsize, uint32_t k)
{
  assert(l>3 && l<=sizeof(HWord)*8 && k>3 && binsize>0 && l<k && JH_NRBITS<32);
  assert(jh_bitmask(HWORDBITS) == HWORD_MAX);

  size_t keylen = (k > l ? k - l : 0) + JH_NRBITS + 1;
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
  fprintf(fh, "  l: %zu, k: %zu, keylen: %zu, nbins: %zu, binsize: %zu, word size: %i\n",
          jh->l, jh->k, jh->keylen, jh->nbins, jh->binsize, (int)HWORDBITS);
  fprintf(fh, "  occupancy: %zu/%zu (%.2f%%) mem: %zu bytes\n",
          jh->nentries, capacity, 100 * (double)jh->nentries / capacity,
          jh->dwords*sizeof(HWord));

  if(jh->nentries > 0) {
    fprintf(fh, "collisions:\n");
    for(j = (1UL<<JH_NRBITS)-1; jh->collisions[j] == 0 && j > 0; j--);
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
  if(o+lenbits < HWORDBITS) firstmsk |= HWORD_MAX << (o+lenbits);

  // char tmp[500];
  // printf(" o: %zu lenbits: %zu firstmask: %s\n", o, lenbits, hword_to_binary(firstmsk, tmp));
  // printf("da%2zu: %s\n", w, hwords_to_binary((HWord*)data+w, 0, 1, tmp));

  neww = entry[0]<<o;
  do {
    oldw = data[w];
    neww = jh_maskmerge(oldw, neww, firstmsk);
  }
  while(!__sync_bool_compare_and_swap(data+w, oldw, neww));

  // printf("da%2zu: %s\n", w, hwords_to_binary((HWord*)data+w, 0, 1, tmp));

  // size_t nwords = jh_nwords(lenbits);
  // printf("givn: %s\n", hwords_to_binary(entry, 0, nwords, tmp));
  // printf("wrte: %s\n", hwords_to_binary((HWord*)data+w, o, nwords, tmp));

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

static inline void _jelly_shuffle(HWord *hash, size_t i, size_t j)
{
  // printf("shuffle %i %i\n", (int)i, (int)j);
  HWord mix;
  hash[i] = jh_mix(hash[i]);
  mix = (hash[i] >> JH_HALFBITS) | (hash[j] << JH_HALFBITS);
  mix = jh_mix(mix);
  hash[i] = (hash[i] & jh_botmsk) | (mix << JH_HALFBITS);
  hash[j] = (hash[j] & jh_topmsk) | (mix >> JH_HALFBITS);
}

static inline void _jelly_unshuffle(HWord *hash, size_t i, size_t j)
{
  // printf("unshuffle %i %i\n", (int)i, (int)j);
  HWord mix;
  mix = (hash[i] >> JH_HALFBITS) | (hash[j] << JH_HALFBITS);
  mix = jh_unmix(mix);
  hash[i] = (hash[i] & jh_botmsk) | (mix << JH_HALFBITS);
  hash[j] = (hash[j] & jh_topmsk) | (mix >> JH_HALFBITS);
  hash[i] = jh_unmix(hash[i]);
}

// klbits is number of bits in top word
static inline void _jelly_hash_key(HWord *hash, size_t kwords, size_t klbits)
{
  size_t i, endbits = JH_HALFBITS + klbits, limit = (klbits < JH_HALFBITS ? 2 : 1);
  HWord mix;

  if(kwords > 1)
  {
    // Mix between words
    for(i = 0; i+limit < kwords; i++)
      _jelly_shuffle(hash, i, i+1);

    if(i > 0)
      _jelly_shuffle(hash, 0, i);

    // Mix last two words
    if(klbits < JH_HALFBITS)
    {
      // printf("partial mix %i\n", (int)i);
      hash[i] = jh_mix(hash[i]);
      mix = (hash[i] >> JH_HALFBITS) | (hash[i+1] << JH_HALFBITS);
      mix = jh_mix_mask(mix, endbits, jh_bitmask(endbits));
      hash[i]   = (hash[i]  & jh_botmsk) | (mix << JH_HALFBITS);
      hash[i+1] = mix >> JH_HALFBITS;
    }
  }

  hash[kwords-1] = jh_mix_mask(hash[kwords-1], klbits, jh_bitmask(klbits));
}

static inline void _jelly_unhash_key(HWord *hash, size_t kwords, size_t klbits)
{
  size_t i, endbits = JH_HALFBITS + klbits;
  HWord mix;

  hash[kwords-1] = jh_unmix_mask(hash[kwords-1], klbits, jh_bitmask(klbits));
  if(kwords == 1) return;

  size_t num_mixes = kwords - 1;

  // unmix last two words
  if(klbits < JH_HALFBITS)
  {
    i = --num_mixes;
    // printf("partial unmix %i\n", (int)i);
    mix = (hash[i] >> JH_HALFBITS) | (hash[i+1] << JH_HALFBITS);
    mix = jh_unmix_mask(mix, endbits, jh_bitmask(endbits));
    hash[i]   = (hash[i] & jh_botmsk) | (mix << JH_HALFBITS);
    hash[i+1] = mix >> JH_HALFBITS;
    hash[i] = jh_unmix(hash[i]);
  }

  if(num_mixes > 0)
    _jelly_unshuffle(hash, 0, num_mixes);

  // unmix between words
  for(i = num_mixes; i > 0; i--)
    _jelly_unshuffle(hash, i-1, i);
}

static inline HKey jelly_hash_find(JellyHash *jh, const char *key,
                                   int insert, int *inserted)
{
  volatile HWord *restrict const data = jh->data;
  const size_t kwords = jh->kwords;
  HWord hash[kwords], entry[jh->hwords+1], found[jh->hwords+1];
  size_t i, w, o, rehashes, binend, pos, lockpos, lockw, locko;
  size_t wordbits, lrem, krem, lbits, kbits;
  HKey loc;
  const HWord emptymsk = jh_bitmask(JH_NRBITS)<<1;
  const HWord lastmsk = jh_bitmask(jh->klbits);

  *inserted = 0;

  // char tmpstr[500];
  // printf("inpt: %s\n", hwords_to_binary(key, 0, kwords, tmpstr));

  // hash
  // for(i = 0; i+1 < kwords; i++) hash[i] = jh_mix(key[i]);
  // hash[i] = jh_mix_mask(key[i] & lastmsk, jh->klbits, lastmsk);
  memcpy(hash, key, kwords*sizeof(HWord));
  hash[kwords-1] &= lastmsk;
  _jelly_hash_key(hash, kwords, jh->klbits);

  // printf("hash: %s\n", hwords_to_binary(hash, 0, kwords, tmpstr));
  // printf("emsk: %s\n", hword_to_binary(emptymsk, tmpstr));

  for(rehashes = 0; rehashes+1 < (1UL<<JH_NRBITS); rehashes++)
  {
    // memset(entry, 0, jh->hwords*sizeof(HWord));
    entry[0] = (rehashes+1) << 1; // collisions id = 1.., lock = 0
    loc = 0;

    // split hash into loc and remaining k-l bits for entry
    HWord word;
    w = 0, o = JH_NRBITS+1;
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

      // printf("searching %zu\n", pos);

      while(1)
      {
        // Wait for lock bit (spin lock)
        while((data[lockw]>>locko) & 1UL);

        _jelly_read(jh, lockw, locko, found, jh->keylen);
        // printf("fnd : %s\n", hword_to_binary(found[0], tmpstr));

        if(memcmp(found, entry, jh->hwords*sizeof(HWord)) == 0) return pos;
        else if((found[0] & emptymsk) * !(found[0] & 1UL)) break;// !empty,!locked,!match
        else {
          // Empty or locked slot - entry not found
          // printf("Found empty!\n");
          if(!insert) return JHASH_NULL;

          if(_jelly_acquire_lock(data, lockw, locko))
          {
            entry[0] |= 1UL; // Set lock = 1

            // Got lock - check still empty
            _jelly_read(jh, lockw, locko, found, jh->keylen);
            if(found[0] & (jh_bitmask(JH_NRBITS)<<1))
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
            __sync_synchronize(); // DEV: is this needed?
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
    // printf("rehash\n");
    // for(i = 0; i+1 < kwords; i++) hash[i] = jh_mix(hash[i]);
    // hash[i] = jh_mix_mask(hash[i], jh->klbits, lastmsk);
    _jelly_hash_key(hash, kwords, jh->klbits);
  }

  return JHASH_NULL;
}


static inline void jelly_hash_get_key(const JellyHash *jh, HKey loc, char *key)
{
  const size_t kwords = jh->kwords;
  size_t bitpos, bitw, bito, i, w, o, nrehashes;
  size_t loccpy, wordbits, lbits, kbits, lrem, krem;
  HWord found[jh->hwords+1], hash[kwords];

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
  o = 1+JH_NRBITS;
  lrem = jh->l; krem = jh->k;

  for(i = 0; i < kwords; i++)
  {
    wordbits = (i+1 == kwords ? jh->klbits : HWORDBITS);
    lbits = (lrem*wordbits+krem-1) / krem;
    kbits = wordbits-lbits;
    // printf("%i kbits:%i lbits:%i\n", (int)i, (int)kbits, (int)lbits);
    hash[i] = (jh_getword(found, w, o) << lbits) | (loccpy & jh_bitmask(lbits));
    loccpy >>= lbits;
    o += kbits; w += (o >= HWORDBITS); o = jh_rembits(o);
    lrem -= lbits; krem -= wordbits;
  }

  // printf("hash : %s\n", hwords_to_binary(hash, 0, kwords, tmp));

  // unhash
  nrehashes = ((found[0]>>1) & jh_bitmask(JH_NRBITS)) - 1;
  // printf("nrehashes: %zu\n", nrehashes);

  // size_t j;
  // HWord mask = jh_bitmask(jh->klbits);

  for(i = 0; i <= nrehashes; i++) {
    // for(j = 0; j+1 < kwords; j++) hash[j] = jh_unmix(hash[j]);
    // hash[j] = jh_unmix_mask(hash[j], jh->klbits, mask);
    _jelly_unhash_key(hash, kwords, jh->klbits);
  }

  // printf("uhsh: %s\n", hwords_to_binary(hash, 0, kwords, tmp));
  memcpy(key, hash, sizeof(HWord)*kwords);
}

#endif /* JELLY_H_ */
