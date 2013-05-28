#ifndef NA63_SIMULATION_PROPAGATE_H
#define NA63_SIMULATION_PROPAGATE_H

#include "Simulation/Track.hh"
#include "Simulation/Simulator.hh"
#include "Simulation/Particle.hh"
#include "Geometry/Material.hh"
#include <curand_kernel.h>

#define TRACK_KEY_DEAD MAX_INT_VALUE - 1
#define TRACK_KEY_AVAILABLE MAX_INT_VALUE

namespace na63 {

typedef struct {
  unsigned N;
  unsigned steps;
  Float dl;
  GPUTrack *tracks;
  int *keys;
  int n_volumes;
  MaterialPars *materials;
  ParticlePars *particles;
  InsideFunction *volume_types;
  VolumePars   *volumes;
  curandState *rng_states;
  unsigned long rng_seed;
} KernelPars;

void PropagateGPU(Simulator *simulator);

} // End namespace na63

#endif /* NA63_SIMULATION_PROPAGATE_H */