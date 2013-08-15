C implementation of the JellyFish hash table
============================================
Isaac Turner  
12 August 2013  
License: Public Domain  

Thread-safe hashtable that uses compressed bit representation.
Fast, low memory and multithreaded. Inspired by the jellyfish kmer-counter
[http://www.cbcb.umd.edu/software/jellyfish].

*Note:* this is experimental code beware changes and bugs.

About
-----

We reduce the memory footprint by using a reversible hash function to allocate
a position in the hash table. For any given element in the hash table, we can
infer information about it by using the reverse hash function on it's position.
This reduces the amount of information (bits) we have to use to store each item.
Additionally, items are stored in a packed bit array.
This implentation is a single .h file and allows multiple threads to add to the
hash at the same time safely.  Should work on 32 and 64 bit systems.

Differences from Jellyfish:
1. single .h file in C instead of C++
2. uses reversible static hash function instead of generating random matrix,
   calculating its inverse and doing matrix multiplication
3. unlimited bound on size of k (number of bits per element) vs 62 bits in JellyFish
4. pt 3 requires us to introduce one bit per element to be used as a write lock

Developed for use with very large hash tables in memory on a single machine with
many cores. This has applications in genome assembly where de Bruijn graphs are
often memory intensive (~100GB of RAM).

Usage
-----

    #include "jellyhash.h"

    ...


    JellyHash jhash;
    jelly_hash_alloc(&jhash, 7, 20, 16); // 20*2^7 entries, 16 bits each

    uint64_t in, out;
    HKey hpos;

    // Find or insert beef
    in = 0xbeef
    hpos = jelly_hash_find(jhash, &in, 1, &inserted);
    if(hpos == HASH_NULL) exit(-1); // hash table full
  
    // Where's the beef? Find and get value from key
    hpos = jelly_hash_find(jhash, &in, 0, &inserted);
    if(hpos == HASH_NULL) exit(-1); // entry not found
    jelly_hash_get_key(jhash, hpos, &out); // get value for `key`
    printf("Found value at position %i: %i\n", (int)key, (int)out);

    jelly_hash_dealloc(&jhash);

To compile copy the two .h files `jellyhash.h` and `twang.h` into your working
directory.


API
---

    void jelly_hash_alloc(JellyHash *jhash, uint32_t l, uint32_t b, uint32_t k)

Allocate a new hash table, with b*2^l entries, bucket size of b, storing values
of size k-bits.  Memory usage is (k-l)*b*2^l. k must be > l, if you want k <= l,
you do not need a hash table and can just use a bit array with no collisions.

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
1. benchmark
2. add support for k < l (just a bit array)
3. faster non-threadsafe functions
4. better 32 bit mix functions
