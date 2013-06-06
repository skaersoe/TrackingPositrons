#ifndef NA63_SIMULATOR
#define NA63_SIMULATOR

#include "types.hh"
#include "box.hh"

namespace na63 {

  void PropagateStraight(particle part);
  void RungeKutta(particle part);
  void MultipleScattering (particle part);
  void EnergyLoss(particle part);

}

#endif