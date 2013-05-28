#include <curand_kernel.h>

#include "Geometry/LibraryCUDA.cuh"
#include "Geometry/Material.hh"
#include "Simulation/Particle.hh"
#include "Simulation/Track.hh"

namespace na63 {

__device__
void CUDA_Bremsstrahlung(GPUTrack& mother, const MaterialPars& material,
    const ParticlePars& particle, const Float dl, curandState *rng_state);

} // End namespace na63