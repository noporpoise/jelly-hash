#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <unistd.h>
#include <pthread.h>
#include <time.h> // srand

// Set number of collisions allowed (max is: 2^NRBITS-1)
// #define NRBITS 5

// Set word size to be used by the hash table
// #define HWord uint32_t
// #define HWORD_MAX UINT32_MAX
// #define HWORDBITS 32

#include "jellyhash.h"

static const char usage[] =
"usage: speedtest [Options] <num_ops>\n"
"  Test hash table speed.  Table size is b*2^l\n"
"    -k <k>  Element size\n"
"    -l <l>  Bits for bucket addressing\n"
"    -b <b>  Number of elements per bucket\n"
"    -t <t>  Number of threads to use\n";

#define die(fmt, ...) call_die(__FILE__, __LINE__, fmt, ##__VA_ARGS__)

void call_die(const char *file, int line, const char *fmt, ...)
__attribute__((format(printf, 3, 4)))
__attribute__((noreturn));

void call_die(const char *file, int line, const char *fmt, ...)
{
  va_list argptr;
  fflush(stdout);
  fprintf(stderr, "[%s:%i] Error: ", file, line);
  va_start(argptr, fmt);
  vfprintf(stderr, fmt, argptr);
  va_end(argptr);
  if(*(fmt + strlen(fmt) - 1) != '\n') fputc('\n', stderr);
  exit(EXIT_FAILURE);
}

// Get a random binary kmer -- useful for testing
static inline void kmer_random(uint32_t kbits, HWord *kmer)
{
  // printf("kbits: %u %u\n", kbits, jh_rembits(kbits));
  size_t i, kwords = jh_nwords(kbits);
  for(i = 0; i < kwords; i++) {
    #if HWORD_MAX > UINT32_MAX
      kmer[i] = (((HWord)rand())<<32) | rand();
    #else
      kmer[i] = rand();
    #endif
  }
  kmer[kwords-1] &= jh_bitmask(kbits);
}

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
  HWord bkmer[kwords];
  int found;
  HKey hpos;

  for(i = 0; i < num_ops; i++)
  {
    // kmer_random(jhash->k, bkmer);
    bkmer[0] = num_ops*thread->threadid + i;
    // printf("kmer: %s\n", hwords_to_binary(bkmer, kwords, tmp));
    hpos = jelly_hash_find(jhash, bkmer, 1, &found);
    if(hpos == HASH_NULL) {
      jelly_hash_print_stats(jhash, stdout);
      die("Hash full!");
    }
    // HWord result[kwords+1];
    // char tmp[kwords*64+1];
    // jelly_hash_get_key(jhash, hpos, result);
    // if(memcmp(result,bkmer,kwords*sizeof(HWord))) {
    //   printf("kmer: %s\n", hwords_to_binary(bkmer, 0, kwords, tmp));
    //   printf("rslt: %s\n", hwords_to_binary(result, 0, kwords, tmp));
    //   die("Mismatch!");
    // }
  }

  pthread_exit(NULL);
}

int main(int argc, char **argv)
{
  uint32_t k = 62, l = 20, b = 32, num_of_threads = 1;
  size_t num_ops;
  int c;

  while ((c = getopt(argc, argv, "k:l:b:t:")) >= 0)
    if (c == 'k') k = atoi(optarg);
    else if (c == 'l') l = atoi(optarg);
    else if (c == 'b') b = atoi(optarg);
    else if (c == 't') num_of_threads = atoi(optarg);

  if(optind == argc) {
    fputs(usage, stderr);
    exit(EXIT_FAILURE);
  }

  num_ops = atol(argv[argc-1]);
  num_ops /= num_of_threads;

  /* initialize random seed: */
  srand(time(NULL) + getpid());

  printf("JellyHash Test k: %u, l: %u, b: %u, threads: %u\n",
         k, l, b, num_of_threads);

  JellyHash jhash;
  jelly_hash_alloc(&jhash, l, b, k);

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

  for(i = 0; i < num_of_threads; i++)
    if(pthread_create(&threads[i], NULL, speedtest, &threaddata[i]) != 0)
      die("Creating thread failed\n");

  for(i = 0; i < num_of_threads; i++)
    if(pthread_join(threads[i], NULL) != 0)
      die("Creating thread failed\n");

  jelly_hash_print_stats(&jhash, stdout);
  jelly_hash_dealloc(&jhash);

  return EXIT_SUCCESS;
}
