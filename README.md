C implementation of the JellyFish hash table
Isaac Turner
12 August 2013
License: Public Domain

Thread-safe hashtable that uses compressed bit representation.
Fast, low memory and multithreaded. Inspired by the jellyfish kmer-counter
[http://www.cbcb.umd.edu/software/jellyfish].

Usage
-----

    #include "jellyhash.h"

    ...


    JellyHash jhash;
    jelly_hash_alloc(&jhash, 7, 16); // 2^7 entries, 16 bits each

    uint64_t in, out;
    HKey key;

    // Find or insert beef
    in = 0xbeef
    key = jelly_hash_find(jhash, &in, 1, &inserted);
    if(key == HASH_NULL) exit(-1); // hash table full
  
    // Where's the beef? Find and get value from key
    key = jelly_hash_find(jhash, &in, 0, &inserted);
    if(key == HASH_NULL) exit(-1); // entry not found
    jelly_hash_get_key(jhash, f, out); // get value for `key`
    printf("Found value at position %i: %i\n", (int)key, (int)out);

    jelly_hash_dealloc(&jhash);


API
---

    void jelly_hash_alloc(JellyHash *jhash, uint32_t l, uint32_t k)

Allocate a new hash table, with 2^l entries, storing values of size k-bits.

    void jelly_hash_dealloc(JellyHash *jhash)

Free the hash table

    HKey jelly_hash_find(JellyHash *jhash, const uint64_t *key,
                         int insert, int *inserted)

Find or insert an item into the hash table.  If an element is not already in the
hash table: returns HASH_NULL if `insert == 0`; otherwise will attempt to
insert the element.  If rehash limit (16) is hit whilst attempting to insert,
returns HASH_NULL.  `key` must be at least `k` bits long. If insertion is
successful, `*inserted` is set to `1`.

    void jelly_hash_get_key(JellyHash *jhash, HKey loc, uint64_t *key)

Use the location of an element to fetch its value.


Cite
----
I am not associated with Guillaume Marçais or Carl Kingsford, the authors
of JellyFish. If you use this software you may wish to cite their paper:

* A fast, lock-free approach for efficient parallel counting of occurrences of k-mers
  Vol. 27 no. 6 2011, pages 764–770, BIOINFORMATICS


TODO
----
1) benchmark
2) add support for k < l (just a bit array)
3) better rehashing?
