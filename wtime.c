#include "wtime.h"

#define _POSIX_C_SOURCE 199309L
#include <time.h>

float wtime(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC_RAW, &ts);

    return 1e-9f * ts.tv_nsec + (float)ts.tv_sec;
}
