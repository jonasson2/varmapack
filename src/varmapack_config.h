#ifndef VARMAPACK_CONFIG_H
#define VARMAPACK_CONFIG_H

#if defined(_WIN32)
  #define HIDDEN
#elif defined(__GNUC__) || defined(__clang__)
  #define HIDDEN __attribute__((visibility("hidden")))
#else
  #define HIDDEN
#endif

#if defined(VARMAPACK_STATIC)
  #define VARMAPACK_API
#elif defined(_WIN32)
  #if defined(VARMAPACK_BUILD)
    #define VARMAPACK_API __declspec(dllexport)
  #else
    #define VARMAPACK_API __declspec(dllimport)
  #endif
#elif defined(__GNUC__) || defined(__clang__)
  #define VARMAPACK_API __attribute__((visibility("default")))
#else
  #define VARMAPACK_API
#endif

#endif
