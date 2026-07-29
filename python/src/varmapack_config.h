#ifndef VARMAPACK_CONFIG_H
#define VARMAPACK_CONFIG_H

#include <stdint.h>
#include <stdbool.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if !defined(__STDC_VERSION__) || __STDC_VERSION__ < 201112L
#error "_Static_assert requires C11"
#endif

_Static_assert(sizeof(int) == 4, "varmapack requires 32-bit int");
_Static_assert(sizeof(void*) == 8, "varmapack requires 64-bit pointers");
_Static_assert(sizeof(long long) == 8, "varmapack requires 64-bit long long");

#if defined(_MSC_VER)
  #define THREADLOCAL __declspec(thread)
#elif defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201112L && \
      !defined(__STDC_NO_THREADS__)
  #define THREADLOCAL _Thread_local
#elif defined(__GNUC__) || defined(__clang__)
  #define THREADLOCAL __thread
#else
  #define THREADLOCAL
#endif

#if defined(_WIN32)
  #define HIDDEN
#elif defined(__GNUC__) || defined(__clang__)
  #define HIDDEN __attribute__((visibility("hidden")))
#else
  #define HIDDEN
#endif

#define TOLOWER(c) (((c) >= 'A' && (c) <= 'Z') ? ((c)-'A'+'a') : (c))
#define STRSET(dst, src) snprintf((dst), sizeof(dst), "%s", (src) ? (src) : "")
#define STRSETF(dst, fmt, ...) snprintf((dst), sizeof(dst), (fmt), __VA_ARGS__)
#define LEN(a) (int)((sizeof(a)/sizeof((a)[0])))
#define CLEAR(dst) memset((dst), 0, sizeof(dst))
#define ALLOC(ptr, count) (((ptr) = calloc((count), sizeof *(ptr))) != 0)
#define FREE(p) do { free(p); (p) = 0; } while (0)

static inline bool intProduct(int a, int b, int *product) {
  if (a < 0 || b < 0 || (a != 0 && b > INT_MAX/a)) return false;
  *product = a*b;
  return true;
}

static inline bool sizeProduct(size_t a, size_t b, size_t *product) {
  if (a != 0 && b > SIZE_MAX/a) return false;
  *product = a*b;
  return true;
}

#endif
