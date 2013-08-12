#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <inttypes.h>
#include <pthread.h> // testing
#include <time.h> // time() for srand()
#include <unistd.h> // getpid() for srand()

#include "jellyhash.h"

#define die(fmt, ...) call_die(__FILE__, __LINE__, fmt, ##__VA_ARGS__)

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

char* uint64s_to_binary(const uint64_t *b, size_t l, char *str)
{
  size_t i, j, k = 0;
  for(i = l-1; i < SIZE_MAX; i--)
    for(j = 0; j < 64; j++)
      str[k++] = '0' + !!(b[i] & (1UL<<(63-j)));
  str[k] = '\0';
  return str;
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
/*
static void test_twang()
{
  uint64_t i, n = 1000000000UL, word, hash, unhash;
  int r = rand();
  int shift = 64 * ((long double)r / (double)RAND_MAX);
  shift=59;
  uint64_t mask = (1UL << shift)-1;
  char tmp[100];
  printf("shift: %i %i %i\nmask: %s\n", shift, r, RAND_MAX, uint64_to_binary(mask, tmp));
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
*/

void* worker(void *ptr)
{
  JellyHash *jhash = (JellyHash*)ptr;
  uint64_t in[10], out[10]; char c; int inserted;
  HKey f;
  char tmp[500];

  sleep(1);
  for(c = 'A'; c <= 'z'; c++) {
    memset(in, c, sizeof(in)); // in[0]=in[1]=c;
    memset(out, 0, sizeof(out));
    f = jelly_hash_find(jhash, in, 1, &inserted);
    if(f == HASH_NULL) die("Hash full");
    jelly_hash_get_key(jhash, f, out);
    printf("%s: %4zu %s\n", (char*)out, (size_t)f, uint64_to_binary(out[0], tmp));
    // printf("%s: %zu\n", out, (size_t)f);
  }
  pthread_exit(NULL);
}

int main()
{
  // uint64_t j, k, ans = 0;
  // for(j = 0; j < 100000000; j++)
  //   for(k = 0; k <= 64; k++)
  //     ans += bitmask(k);
  // printf("k: %zu j: %zu ans: %zu\n", (size_t)k, (size_t)j, (size_t)ans);
  // exit(-1);

  // Initiate rand
  srand(time(NULL) + getpid());
  // test_twang();

  JellyHash jhash;
  jelly_hash_alloc(&jhash, 7, 8*9);

  // worker((void*)&jhash);

  int rc, i, num_of_threads = 1000;
  pthread_t threads[num_of_threads];

  for(i = 0; i < num_of_threads; i++)
  {
    rc = pthread_create(&threads[i], NULL, worker, (void*)&jhash);
    if(rc != 0) die("Creating thread failed");
  }

  for(i = 0; i < num_of_threads; i++) {
    rc = pthread_join(threads[i], NULL);
    if(rc != 0) die("Joining thread failed");
  }

  jelly_hash_dealloc(&jhash);

  return EXIT_SUCCESS;
}
