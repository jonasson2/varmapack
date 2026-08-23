// TimeUtil.h: shared timing utilities for C benchmarks.

#ifndef TIMEUTIL_H
#define TIMEUTIL_H

#include <stdbool.h>
#include <stdint.h>

typedef struct {
  uint64_t start;
  uint64_t deadline;
  uint64_t now;
  int reps;
} time_loop;

uint64_t clock_nsec(void);
void consume_double(double x);
void warmup_cpu(double seconds);
void time_loop_start(time_loop *timer, double seconds);
bool time_loop_next(time_loop *timer);
double time_loop_average(time_loop *timer, uint64_t elapsed);
double time_loop_nsec_per(time_loop *timer, int values);

#endif
