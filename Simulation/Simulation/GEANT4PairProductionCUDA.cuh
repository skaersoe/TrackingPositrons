#ifndef NA63_SIMULATION_GEANT4PAIRPRODUCTION_CUH
#define NA63_SIMULATION_GEANT4PAIRPRODUCTION_CUH

#include "Geometry/LibraryCUDA.cuh"

namespace na63 {

__device__
void CUDA_GEANT4PairProduction(
    GPUTrack* mother,
    const ParticlePars* particle,
    const MaterialPars* material,
    const Float dl,
    curandState *rng_state, 
    const int index
);

}

#endif /* NA63_SIMULATION_GEANT4PAIRPRODUCTION_CUH */