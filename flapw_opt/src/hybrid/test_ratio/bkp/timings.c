#include "timings.h"


int read_clock(struct timeval *t)
{
    return gettimeofday(t, NULL);
}

int elapsed_time(struct timeval *start, struct timeval *end)
{
    return (int)(end->tv_sec  - start->tv_sec) * 1e6 + 
		   (int)(end->tv_usec - start->tv_usec);
}