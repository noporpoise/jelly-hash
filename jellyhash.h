#ifndef JELLY_H_
#define JELLY_H_

#include <inttypes.h>

// rehash limit is 2^RBITS-1 -> rbits=4 => 15
#define RBITS 4

typedef struct
{
  volatile uint64_t *const data;
  const uint32_t l, k, keylen, hwords; // hwords is num of words of k in hash
  const uint32_t kwords, klbits; // klbits is bits in top word
  const uint32_t strt_lbits, full_lbits, last_lbits;
  const uint32_t strt_kbits, full_kbits, last_kbits;
  const uint64_t nkeys, lmask;
} JellyHash;

typedef uint64_t HKey;
#define HASH_NULL UINT64_MAX

// size is 2^l, k is kmer-size
void jelly_hash_alloc(JellyHash *jhash, uint32_t l, uint32_t k);
void jelly_hash_dealloc(JellyHash *jhash);

HKey jelly_hash_find(JellyHash *jhash, uint64_t *key, int insert, int *inserted);
// HKey jelly_hash_find(JellyHash *jhash, uint64_t *key);
// HKey jelly_hash_find_or_insert(JellyHash *jhash, uint64_t *key, int *found);
void jelly_hash_get_key(JellyHash *jhash, HKey pos, uint64_t *key);

#endif /* JELLY_H_ */
