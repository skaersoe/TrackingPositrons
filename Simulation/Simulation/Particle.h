#ifndef NA63_SIMULATION_PARTICLE_H
#define NA63_SIMULATION_PARTICLE_H

namespace Simulation {

  typedef struct {
    float r[3];
    float p[3];
    int id;
    // Align to 32 bytes
    char padding[4];
  } simple_particle_t;

}

#endif