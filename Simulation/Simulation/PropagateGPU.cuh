#ifndef NA63_SIMULATION_PROPAGATE_H
#define NA63_SIMULATION_PROPAGATE_H

#include "Simulation/Track.hh"
#include "Simulation/Simulator.hh"
#include "Simulation/Particle.hh"
#include "Geometry/Material.hh"
#include <curand_kernel.h>

#define TRACK_KEY_FREE MAX_INT_VALUE
#define TRACK_KEY_DEAD TRACK_KEY_FREE - 1
#define TRACK_KEY_WAITING TRACK_KEY_FREE - 2

namespace na63 {

typedef struct {
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
  SortMethod sorting;
} KernelPars;

void PropagateGPU(Simulator *simulator);

} // End namespace na63

#endif /* NA63_SIMULATION_PROPAGATE_H */