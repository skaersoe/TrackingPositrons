#ifndef NA63_SIMULATION_PROPAGATE_H
#define NA63_SIMULATION_PROPAGATE_H

#include "Simulation/Track.h"
#include "Simulation/Arguments.h"

typedef struct {
  unsigned N;
  unsigned steps;
  float dt;
} KernelPars;

namespace na63 {

  void Propagate(Track *t, SimulatorPars args);

}

#endif