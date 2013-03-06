#ifndef PROPAGATE_H
#define PROPAGATE_H

#include "Simulation/Particle.h"
#include "Simulation/Arguments.h"

namespace Simulation {

  void propagate(simple_particle_t *p, simulator_args_t args);

}

#endif