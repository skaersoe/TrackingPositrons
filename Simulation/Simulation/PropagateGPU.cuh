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
    Float dl;
    GPUTrack *tracks;
    int *keys;
    MaterialPars *materials;
    ParticlePars *particles;
    VolumePars   *volumes;
  } KernelPars;

  void PropagateGPU(Simulator *simulator);

}

#endif /* NA63_SIMULATION_PROPAGATE_H */