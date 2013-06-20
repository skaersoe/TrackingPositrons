#include "Simulation/BetheEnergyLossCUDA.cuh"
#include "Simulation/TrackGPU.cuh"
#include "Simulation/LandauCUDA.cuh"
#include "Geometry/LibraryCUDA.cuh"

// K/A = 0.307075 MeV g^-1 cm^2

namespace na63 {

__device__
void CUDA_BetheEnergyLoss(GPUTrack& track, const ParticlePars& particle,
    const MaterialPars& material, const Float dl, const Float index,
    curandState *rng_state) {

  // Only treat living particles with charge
  if (track.state != ALIVE || track.charge == 0) return;

  Float mass = particle.mass;

  // Don't treat massless particles
  if (mass == 0) return;

  // Get -<dE/dx> and sigma
  LandauParameters p = CUDA_LandauEnergyLossParameters(
      CUDA_Gamma(track.momentum[3],mass),mass,
      material.atomic_number,material.density,
      material.atomic_weight,material.mean_excitation_potential,dl);

  // Get random number from Landau distribution
  Float energy_loss = ThrowLandau(p.mpv,4*p.xi,curand_uniform(rng_state));

  // Update track
  // if (p.mpv != p.mpv) {
  //   printf("CUDA_BetheEnergyLoss sending NaN\n");
  //   return;
  // }

  // if (energy_loss != energy_loss) {
  //   printf("CUDA_BetheEnergyLoss sending NaN\n");
  //   return;
  // }

  CUDA_UpdateEnergy(track.momentum,mass,-energy_loss,index);
}

} // End namespace na63