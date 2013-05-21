#include "Simulation/BetheEnergyLoss.hh"
#include "Geometry/Library.hh"

// K/A = 0.307075 MeV g^-1 cm^2

namespace na63 {

void BetheEnergyLoss(Track& track, const Material& material,
    const Float dl) {

  TRandom3 rng;

  rng.SetSeed((size_t)&track+(size_t)BetheEnergyLoss);

  // Only treat particles with charge
  if (track.charge() == 0) return;

  // Get -<dE/dx> and sigma
  LandauParameters p = GetSkewedLandauParameters(
      track.beta(),track.mass(),track.charge(),material.atomic_number(),
      material.mean_excitation_potential(),dl);

  // Get random number from Landau distribution
  Float energy_loss = rng.Landau(p.mean,p.sigma);

  // Update track
  track.UpdateMomentum(-energy_loss * dl);
}

} // End namespace na63