#include "Simulation/GetTime.hh"

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

double InSeconds(const timespec ts) {
  double s(ts.tv_sec);
  return s + (1e-9*(double)ts.tv_nsec);
}

} // End namespace na63