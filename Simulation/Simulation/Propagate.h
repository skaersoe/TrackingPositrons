#ifndef NA63_SIMULATION_PROPAGATE_H
#define NA63_SIMULATION_PROPAGATE_H

#include "Simulation/Particle.h"
#include "Simulation/Arguments.h"

typedef struct {
  unsigned N;
  unsigned steps;
  float dt;
} kernel_args_t;

namespace Simulation {

  void propagate(simple_particle_t *p, simulator_args_t args);

}

#endif