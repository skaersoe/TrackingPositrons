#ifndef GENERATEPARTICLES_H
#define GENERATEPARTICLES_H

// #include "HepMC/GenParticle.h"
#include <cuda.h>
#include <cuda_runtime.h>

#include "Simulation/Simulator.h"
#include "Simulation/Particle.h"
#include "Simulation/Arguments.h"

namespace Simulation {

  void generateParticles(simple_particle_t* p, simulator_args_t args);

}

#endif