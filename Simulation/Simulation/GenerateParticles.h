#ifndef GENERATEPARTICLES_H
#define GENERATEPARTICLES_H

// #include "HepMC/GenParticle.h"
#include <cuda.h>
#include <cuda_runtime.h>

#include "Particle.h"

void GenerateParticles(simple_particle_t* p, const int N);

#endif