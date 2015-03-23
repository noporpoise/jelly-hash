#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <assert.h>

/**
 * Return number of digits required to represent `num` in base 10.
 * Examples:
 *   num_of_digits(0)   = 1
 *   num_of_digits(1)   = 1
 *   num_of_digits(10)  = 2
 *   num_of_digits(123) = 3
 */
size_t num_of_digits(size_t num)
{
  size_t digits = 1;
  while(1) {
    if(num < 10) return digits;
    if(num < 100) return digits+1;
    if(num < 1000) return digits+2;
    if(num < 10000) return digits+3;
    num /= 10000;
    digits += 4;
  }
  return digits;
}

// result must be long enough for result + 1 ('\0'). Max length required is:
// strlen('18,446,744,073,709,551,615')+1 = 27
// returns pointer to result
char* ulong_to_str(unsigned long num, char *result)
{
  unsigned int digits = num_of_digits(num);
  unsigned int i, num_commas = (digits-1) / 3;
  char *p = result + digits + num_commas;
  *(p--) = '\0';

  for(i = 0; i < digits; i++, num /= 10) {
    if(i > 0 && i % 3 == 0) *(p--) = ',';
    *(p--) = '0' + (num % 10);
  }

  return result;
}

// result must be long enough for result + 1 ('\0'). Max length required is:
// strlen('-9,223,372,036,854,775,808')+1 = 27
char* long_to_str(long num, char *result)
{
  if(num < 0) {
    result[0] = '-';
    result++;
    num = -num;
  }

  ulong_to_str((unsigned long)num, result);

  return result;
}

// result must be long enough for result + 1 ('\0').
// Max length required is: 26+1+decimals+1 = 28+decimals bytes
// strlen('-9,223,372,036,854,775,808') = 27
// strlen('.') = 1
// +1 for \0
char* double_to_str(double num, int decimals, char* str)
{
  if(isnan(num)) return strcpy(str, "NaN");
  else if(isinf(num)) return strcpy(str, "Inf");

  unsigned long whole_units = (unsigned long)num;
  num -= whole_units;

  char decstr[2+decimals+1];
  sprintf(decstr, "%.*lf", decimals, num);
  if(decstr[0] == '1') whole_units++;

  ulong_to_str(whole_units, str);

  if(decimals > 0)
  {
    size_t offset = strlen(str);
    strcpy(str+offset, decstr+1);
  }

  return str;
}

// Format a number
static inline char* units_to_str(double num, int decimals, char* str,
                                 const char **units, size_t nunits, size_t usize)
{
  assert(nunits > 0 && usize > 0);

  if(isnan(num)) { sprintf(str, "NaN%s", units[0]); return str; }
  else if(isinf(num)) { sprintf(str, "Inf%s", units[0]); return str; }

  size_t unit = 0;
  double num_tmp = num, num_of_units;

  while(num_tmp >= usize && unit+1 < nunits) { unit++; num_tmp /= usize; }

  num_of_units = num / pow(usize, unit);
  double_to_str(num_of_units, decimals, str);

  char *ptr = str+strlen(str)-1;
  if(decimals > 0) {
    // Trim excess zeros
    while(ptr > str && *ptr == '0') ptr--;
    if(*ptr == '.') ptr--;
  }
  strcpy(ptr+1, units[unit]);

  return str;
}

// str must be 26 + 3 + 1 + num decimals + 1 = 31+decimals bytes
// breakdown:
// strlen('18,446,744,073,709,551,615') = 26
// strlen(' GB') = 3
// strlen('.') = 1
// +1 for '\0'
char* bytes_to_str(unsigned long num, int decimals, char* str)
{
  const char *units[7] = {"B", "KB", "MB", "GB", "TB", "PB", "EB"};
  return units_to_str(num, decimals, str, units, 7, 1024);
}

// Number to string using G to mean 10^9, M to mean 10^6 etc
char* num_to_str(double num, int decimals, char* str)
{
  const char *units[4] = {"", "K", "M", "G"};
  return units_to_str(num, decimals, str, units, 4, 1000);
}
