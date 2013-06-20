#ifndef NA63_SIMULATION_BETHEENERGYLOSSCUDA_CUH
#define NA63_SIMULATION_BETHEENERGYLOSSCUDA_CUH

#include <curand_kernel.h>

#include "Simulation/BetheEnergyLoss.hh"
#include "Simulation/TrackGPU.cuh"
#include "Geometry/LibraryCUDA.cuh"

// K/A = 0.307075 MeV g^-1 cm^2

namespace na63 {

__device__
inline LandauParameters CUDA_LandauEnergyLossParameters(const Float gamma,
    const Float mass, const Float atomic_number,
    const Float density, const Float atomic_weight,
    const Float mean_excitation_potential, const Float dl) {

  LandauParameters p;

  // Calculate necessary values
  Float beta = CUDA_Beta(gamma);
  Float beta_squared = beta*beta;
  Float gamma_squared = gamma*gamma;

  p.xi = 0.5 * dl * 0.1535 * density * atomic_number /
      (atomic_weight * beta_squared);

  p.mpv = p.xi * (std::log(2.0 * mass * beta_squared
      * gamma_squared / mean_excitation_potential)
      + std::log(p.xi/mean_excitation_potential) + 0.200 - beta_squared);

  return p;

}

__device__
inline Float CUDA_Bethe_dEdx(const Float beta,
    const Float mass, const Float charge, const Float atomic_number,
    const Float mean_excitation_potential, const Float density,
    const Float atomic_weight) {

  Float beta_squared = beta*beta;
  Float betagamma_squared = std::pow(Gamma(beta),2) * beta_squared;

  Float mean;

  Float s = kElectronMass / mass;
  Float W_max = 2.0 * kElectronMass * betagamma_squared
      / (1.0 + 2.0 * s * sqrt(1 + betagamma_squared) + s*s);

  mean = 0.1535 * density * atomic_number / atomic_weight * charge*charge
      / beta_squared;

  mean *= std::log((2.0 * kElectronMass * betagamma_squared * W_max)
      / (mean_excitation_potential*mean_excitation_potential))
      - 2.0 * beta_squared;

  return mean;

}

__device__
inline Float CUDA_BetheElectron_dEdx(const Float density,
    const Float atomic_number, const Float beta, const Float kinetic_energy,
    const Float mean_excitation_potential, const Float charge) {

  // Kinetic energy in terms of electron masses
  Float tau = kinetic_energy / kElectronMass;
  Float beta_squared = beta*beta;

  Float F_tau_electron = 1.0 - beta_squared
      + (tau*tau / 8.0 - std::log(2.0) * (2.0*kElectronRadius + 1.0))
      / pow(tau + 1.0,2);

  Float F_tau_positron = 2 * std::log(2) - beta_squared / 12
      * (23 + 14/(tau + 2) + 10/std::pow(tau + 2,2) + 4/std::pow(tau+2,3));

  Float mean = 2.0 * kPi * kAvogadro * kElectronRadiusSquared * kElectronMass
      * density * atomic_number / (beta*beta);
  mean *= std::log((tau*tau * (tau + 2.0))
      / (2.0 * std::pow(1e-6 * mean_excitation_potential / kElectronMass,2)))
      + ((charge < 0) ? F_tau_electron : F_tau_positron);

  return mean;

}

// __device__ inline
// LandauParameters CUDA_GetSkewedLandauParameters(const Float beta,
//     const Float mass, const Float atomic_number,
//     const Float mean_excitation_potential, const Float dl) {

//   LandauParameters p;

//   // Calculate necessary values
//   const Float beta_squared = pow(beta,2);
//   const Float gamma = CUDA_Gamma(beta);
//   const Float gamma_squared = pow(gamma,2);
//   const Float xi = 4.0 * 0.307075 * atomic_number * dl / beta_squared;

//   p.mpv = xi * (log(2 * mass * /*pow(kC,2) **/ beta_squared
//       * gamma_squared / mean_excitation_potential)
//       + log(xi/mean_excitation_potential) + 0.200 - beta_squared);
//   p.xi = 4 * xi;

//   return p;

// }

__device__
void CUDA_BetheEnergyLoss(GPUTrack& track, const ParticlePars& particle,
    const MaterialPars& material, const Float dl, const Float index,
    curandState *rng_state);

} // End namespace na63

#endif /* NA63_SIMULATION_BETHEENERGYLOSSCUDA_CUH */