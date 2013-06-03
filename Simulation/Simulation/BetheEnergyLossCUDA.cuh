#ifndef NA63_SIMULATION_BETHEENERGYLOSSCUDA_CUH
#define NA63_SIMULATION_BETHEENERGYLOSSCUDA_CUH

#include <curand_kernel.h>

#include "Simulation/BetheEnergyLoss.hh"
#include "Simulation/TrackGPU.cuh"
#include "Geometry/LibraryCUDA.cuh"

// K/A = 0.307075 MeV g^-1 cm^2

namespace na63 {

__device__ inline
LandauParameters CUDA_GetSkewedLandauParameters(const Float beta,
    const Float mass, const Float atomic_number,
    const Float mean_excitation_potential, const Float dl) {

  LandauParameters p;

  // Calculate necessary values
  const Float beta_squared = pow(beta,2);
  const Float gamma = CUDA_Gamma(beta);
  const Float gamma_squared = pow(gamma,2);
  const Float xi = 4.0 * 0.307075 * atomic_number * dl / beta_squared;

  p.mean = xi * (log(2 * mass * /*pow(kC,2) **/ beta_squared
      * gamma_squared / mean_excitation_potential)
      + log(xi/mean_excitation_potential) + 0.200 - beta_squared);
  p.sigma = 4 * xi;

  return p;

}

__device__
void CUDA_BetheEnergyLoss(GPUTrack& track, const ParticlePars& particle,
    const MaterialPars& material, const Float dl, curandState *rng_state);

} // End namespace na63

#endif /* NA63_SIMULATION_BETHEENERGYLOSSCUDA_CUH */