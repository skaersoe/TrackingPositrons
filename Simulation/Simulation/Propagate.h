#ifndef NA63_SIMULATION_PROPAGATE_H
#define NA63_SIMULATION_PROPAGATE_H

#include "Simulation/Track.h"
#include "Simulation/Simulator.h"
#include "Simulation/Particle.h"
#include "Geometry/Material.h"

namespace na63 {

  typedef struct {
    unsigned N;
    unsigned steps;
    float dt;
    MaterialPars *material_arr;
    ParticlePars *particle_arr;
  } KernelPars;

  void Propagate(Track *t, SimulatorPars args);

}

#endif