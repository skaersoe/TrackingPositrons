#ifndef STEP_H
#define STEP_H

#include "Particle.h"

void timestep(simple_particle_t* p, float dt);
void timestep(float* r, float* p, float m, float q, float dt);

#endif