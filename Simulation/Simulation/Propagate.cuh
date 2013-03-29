#ifndef NA63_SIMULATION_PROPAGATE_H
#define NA63_SIMULATION_PROPAGATE_H

#include "Simulation/Track.hh"
#include "Simulation/Simulator.hh"
#include "Simulation/Particle.hh"
#include "Geometry/Material.hh"

namespace na63 {

  typedef struct {
    unsigned N;
    unsigned steps;
    float dt;
    MaterialPars *material_arr;
    ParticlePars *particle_arr;
    VolumePars   *volume_arr;
  } KernelPars;

  void Propagate(Track *t, SimulatorPars args);

}

#endif