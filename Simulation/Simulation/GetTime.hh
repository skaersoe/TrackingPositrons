#ifndef NA63_COMMON_GETTIME_H
#define NA63_COMMON_GETTIME_H

// Stolen from http://stackoverflow.com/questions/5167269/clock-gettime-alternative-in-mac-os-x

#include <time.h>
#include <sys/time.h>

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

namespace na63 {

timespec GetTime();

double InSeconds(const timespec ts);

} // End namespace na63

#endif /* NA63_COMMON_GETTIME_H */