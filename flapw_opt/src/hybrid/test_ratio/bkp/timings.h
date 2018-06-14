/*
 * timings.h
 *
 *
 */

#ifndef TIMINGS_H_
#define TIMINGS_H_

#include <sys/time.h>
#include <time.h>


/* Read current time */
int read_clock(struct timeval *t);

/* Calculate duration */
int elapsed_time(struct timeval *start, struct timeval *end);

#endif