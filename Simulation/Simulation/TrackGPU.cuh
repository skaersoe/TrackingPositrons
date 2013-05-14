#ifndef NA63_SIMULATION_TRACKGPU_H
#define NA63_SIMULATION_TRACKGPU_H

#include "Geometry/LibraryCUDA.cuh"
#include "Simulation/Track.hh"
#include "Simulation/Particle.hh"

namespace na63 {

__device__
int VolumeQuery(const GPUFourVector& position, const VolumePars* volumes,
    const InsideFunction* functions, const int volume_index);

__device__
void Step(GPUTrack& track, const ParticlePars& particle, const Float dl);

} // End namespace na63

#endif /* NA63_SIMULATION_TRACKGPU_H */