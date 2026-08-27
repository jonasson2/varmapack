// TimeUtil.c: shared timing utilities for C benchmarks.

#include <math.h>
#include <time.h>
#ifdef _WIN32
#include <windows.h>
#endif
#include "TimeUtil.h"

uint64_t clock_nsec(void) {
#ifdef _WIN32
  LARGE_INTEGER freq, counter;
  QueryPerformanceFrequency(&freq);
  QueryPerformanceCounter(&counter);
  return (uint64_t)(1000000000*(double)counter.QuadPart/(double)freq.QuadPart);
#else
#ifdef CLOCK_MONOTONIC
  struct timespec ts;
  if (clock_gettime(CLOCK_MONOTONIC, &ts) == 0) {
    return (uint64_t)ts.tv_sec*1000000000ull + (uint64_t)ts.tv_nsec;
  }
#endif
  return (uint64_t)(1000000000*(double)clock()/(double)CLOCKS_PER_SEC);
#endif
}

void consume_double(double x) {
  static volatile double sink;
  sink += x;
  (void)sink;
}

void warmup_cpu(double seconds) {
  volatile double sink = 0;
  double x = 1.000001;
  uint64_t t, deadline;
  if (seconds <= 0) return;
  t = clock_nsec();
  deadline = t + (uint64_t)(seconds*1e9);
  while (t < deadline) {
    for (int i=0; i<1000; i++) {
      sink += log(x);
      x += .000001;
      if (x >= 2) x = 1.000001;
    }
    t = clock_nsec();
  }
  (void)sink;
}

void time_loop_start(time_loop *timer, double seconds) {
  timer->start = clock_nsec();
  timer->deadline = timer->start + (uint64_t)(seconds*1e9);
  timer->now = timer->start;
  timer->reps = 0;
}

bool time_loop_next(time_loop *timer) {
  if (timer->reps > 0) timer->now = clock_nsec();
  if (timer->now >= timer->deadline) return false;
  timer->reps++;
  return true;
}

double time_loop_average(time_loop *timer, uint64_t elapsed) {
  return (double)elapsed/timer->reps;
}

double time_loop_nsec_per(time_loop *timer, int values) {
  return (double)(timer->now - timer->start)/((double)timer->reps*values);
}
