#ifndef TIMING_H
#define TIMING_H

#include <sys/time.h>

int read_clock(struct timeval *t);
int elapsed_time(struct timeval *start, struct timeval *end);

#endif
