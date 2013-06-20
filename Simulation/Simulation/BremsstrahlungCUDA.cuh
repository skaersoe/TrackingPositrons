#ifndef NA63_SIMULATION_CUDABREMSSTRAHLUNG_CUH
#define NA63_SIMULATION_CUDABREMSSTRAHLUNG_CUH

#include <curand_kernel.h>

#include "Geometry/LibraryCUDA.cuh"
#include "Geometry/Material.hh"
#include "Simulation/Particle.hh"
#include "Simulation/Track.hh"

namespace na63 {

__device__
void CUDA_Bremsstrahlung(
    GPUTrack& mother,
    const ParticlePars& particle,
    const MaterialPars& material,
    const Float dl,
    curandState *rng_state, 
    const int index);

} // End namespace na63

#endif /* NA63_SIMULATION_CUDABREMSSTRAHLUNG_CUH */