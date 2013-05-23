#ifndef NA63_SIMULATION_TRACKGPU_H
#define NA63_SIMULATION_TRACKGPU_H

#include "Geometry/LibraryCUDA.cuh"
#include "Simulation/Track.hh"
#include "Simulation/Particle.hh"

namespace na63 {

__device__
int VolumeQuery(const GPUTrack& track, const VolumePars* volumes,
    const unsigned n_volumes);

__device__
void Step(GPUTrack& track, const ParticlePars& particle, const Float dl);

__device__
void CUDA_UpdateMomentum(GPUFourVector& momentum, const Float change);

} // End namespace na63

#endif /* NA63_SIMULATION_TRACKGPU_H */