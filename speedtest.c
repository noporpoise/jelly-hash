#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <unistd.h>
#include <pthread.h>
#include <assert.h>
#include <time.h> // srand

#include "num_format.h"

// Set number of collisions allowed (max is: 2^JH_NRBITS-1)
// #define JH_NRBITS 5

// Set word size to be used by the hash table
// #define JH_FORCE_64 1
// #define JH_FORCE_32 1

#include "jellyhash.h"

void print_usage()
{
  printf("usage: speedtest [Options] <num_ops>\n"
"  Test hash table speed.  Table capacity is: b*2^l. Memory is: (k-l+%i)*capacity.\n"
"  Rehash limit is: %i\n"
"    -k <k>  Element size\n"
"    -l <l>  Bits for bucket addressing\n"
"    -b <b>  Number of elements per bucket\n"
"    -t <t>  Number of threads to use\n", JH_NRBITS+1, (1<<JH_NRBITS)-1);
  exit(EXIT_FAILURE);
}

#define die(fmt,...) do { \
  fprintf(stderr, "[%s:%i] Error: %s() "fmt"\n", __FILE__, __LINE__, __func__, ##__VA_ARGS__); \
  exit(EXIT_FAILURE); \
} while(0)


typedef struct {
  int threadid;
  JellyHash *jhash;
  size_t num_ops;
} TestThread;

static inline void* speedtest(void *ptr)
{
  TestThread *thread = (TestThread*)ptr;

  srand(time(NULL) + getpid() + thread->threadid);
  JellyHash *jhash = thread->jhash;
  size_t num_ops = thread->num_ops;

  size_t i, kwords = jh_nwords(jhash->k);
  HWord bkmer[kwords], result[kwords];
  int found;
  HKey hpos;
  memset(bkmer, 0, sizeof(HWord)*kwords);
  memset(result, 0, sizeof(HWord)*kwords);

  for(i = 0; i < num_ops; i++)
  {
    bkmer[0] = num_ops*thread->threadid + i;
    hpos = jelly_hash_find(jhash, (char*)bkmer, 1, &found);
    if(hpos == JHASH_NULL) {
      jelly_hash_print_stats(jhash, stdout);
      die("Hash full!");
    }
  }
  pthread_exit(NULL);
}

int main(int argc, char **argv)
{
  uint32_t k = 62, l = 20, b = 32, num_of_threads = 1;
  size_t num_ops;
  int c;

  while ((c = getopt(argc, argv, "k:l:b:t:h")) >= 0)
    if (c == 'k') k = atoi(optarg);
    else if (c == 'l') l = atoi(optarg);
    else if (c == 'b') b = atoi(optarg);
    else if (c == 't') num_of_threads = atoi(optarg);
    else if (c == 'h') print_usage();

  if(optind == argc)
    print_usage();

  num_ops = atol(argv[argc-1]);
  num_ops /= num_of_threads;

  /* initialize random seed: */
  srand(time(NULL) + getpid());

  JellyHash jhash;
  jelly_hash_alloc(&jhash, l, b, k);

  char memstr[50];
  bytes_to_str((jhash.nbins * jhash.binsize * jhash.keylen) / 8, 1, memstr);
  printf("JellyHash Test k: %u, l: %u, b: %u, threads: %u mem: %s\n",
         k, l, b, num_of_threads, memstr);

  jelly_hash_print_stats(&jhash, stdout);

  // TestThread tt = {.threadid = 101, .jhash = &jhash};
  // speedtest(&tt);

  uint32_t i;
  pthread_t threads[num_of_threads];
  TestThread threaddata[num_of_threads];

  for(i = 0; i < num_of_threads; i++) {
    threaddata[i].threadid = i;
    threaddata[i].jhash = &jhash;
    threaddata[i].num_ops = num_ops;
  }

  if(num_of_threads <= 1) {
    speedtest(&threaddata[0]);
  }
  else {
    for(i = 0; i < num_of_threads; i++)
      if(pthread_create(&threads[i], NULL, speedtest, &threaddata[i]) != 0)
        die("Creating thread failed\n");

    for(i = 0; i < num_of_threads; i++)
      if(pthread_join(threads[i], NULL) != 0)
        die("Creating thread failed\n");
  }

  jelly_hash_print_stats(&jhash, stdout);
  jelly_hash_dealloc(&jhash);

  return EXIT_SUCCESS;
}
