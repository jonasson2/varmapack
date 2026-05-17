#ifndef VARMAPACK_CONFIG_H
#define VARMAPACK_CONFIG_H

#if defined(_WIN32)
  #define HIDDEN
#elif defined(__GNUC__) || defined(__clang__)
  #define HIDDEN __attribute__((visibility("hidden")))
#else
  #define HIDDEN
#endif

#endif
