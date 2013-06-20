#ifndef NA63_SIMULATION_ANGULARDISTRIBUTIONCUDA_CUH
#define NA63_SIMULATION_ANGULARDISTRIBUTIONCUDA_CUH

#include "Geometry/LibraryCUDA.cuh"

namespace na63 {

__device__
void CUDA_ModifiedTsai_SampleDirection(GPUThreeVector& local_direction,
    const GPUTrack* track, const Float mass, curandState *rng_state);

} // End namespace na63

#endif /* NA63_SIMULATION_ANGULARDISTRIBUTIONCUDA_CUH */