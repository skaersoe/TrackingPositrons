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

  Float mass = particle.mass;

  // Don't handle electrons for now
  if (mass < 1 * MeV) return;

  // Get -<dE/dx> and sigma
  LandauParameters p = CUDA_GetSkewedLandauParameters(
      CUDA_Beta(track.momentum[3],mass),mass,
      material.atomic_number,material.mean_excitation_potential,dl);

  // Get random number from Landau distribution
  Float energy_loss = ThrowLandau(p.mean,p.sigma,curand_uniform(rng_state));

  // Update track
  CUDA_UpdateEnergy(track.momentum,mass,-energy_loss * dl);
}

} // End namespace na63