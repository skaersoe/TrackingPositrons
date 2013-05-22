#include "Simulation/BetheEnergyLossCUDA.cuh"
#include "Simulation/TrackGPU.cuh"
#include "Simulation/Landau.cuh"
#include "Geometry/LibraryCUDA.cuh"

// K/A = 0.307075 MeV g^-1 cm^2

namespace na63 {

__device__
void CUDA_BetheEnergyLoss(GPUTrack& track, const ParticlePars& particle,
    const MaterialPars& material, const Float dl, curandState *rng_state) {

  // Only treat particles with charge
  if (track.charge == 0) return;

  // Get -<dE/dx> and sigma
  const LandauParameters p = CUDA_GetSkewedLandauParameters(
      CUDA_Beta(track.momentum[3],particle.mass),particle.mass,track.charge,
      material.atomic_number,material.mean_excitation_potential,dl);

  // Get random number from Landau distribution
  const Float energy_loss = ThrowLandau(p.mean,p.sigma,curand_uniform(rng_state));

  // Update track
  CUDA_UpdateMomentum(track.momentum,-energy_loss * dl);
}

} // End namespace na63