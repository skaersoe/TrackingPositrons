#ifndef GENERATEPARTICLES_H
#define GENERATEPARTICLES_H

// #include "HepMC/GenParticle.h"
#include <cuda.h>
#include <cuda_runtime.h>

#include "Simulation/Particle.h"
#include "Simulation/LaunchArguments.h"

void GenerateParticles(simple_particle_t* p, launch_args_t args);

#endif