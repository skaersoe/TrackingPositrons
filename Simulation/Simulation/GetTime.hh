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

timespec GetTime() {

  struct timespec ts;

  #ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
  clock_serv_t cclock;
  mach_timespec_t mts;
  host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
  clock_get_time(cclock, &mts);
  mach_port_deallocate(mach_task_self(), cclock);
  ts.tv_sec = mts.tv_sec;
  ts.tv_nsec = mts.tv_nsec;

  #else
  clock_gettime(CLOCK_REALTIME, &ts);
  #endif

  return ts;

}

double InSeconds(timespec ts) {
  double s(ts.tv_sec);
  return s + (1e-9*(double)ts.tv_nsec);
}

} // End namespace na63

#endif /* NA63_COMMON_GETTIME_H */