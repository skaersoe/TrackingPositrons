#include "Simulation/Landau.hh"
#include "Simulation/BetheEnergyLoss.hh"
#include "Geometry/Library.hh"

// K/A = 0.307075 MeV g^-1 cm^2

namespace na63 {

TRandom3 rng;

void BetheEnergyLoss(Track& track, const Material& material,
    const Float dl) {

  // Only treat particles with charge
  if (track.charge() == 0) return;

  // Get -<dE/dx> and sigma
  LandauParameters p = LandauEnergyLossParameters(
      track.beta(),track.mass(),material.atomic_number(),
      material.mean_excitation_potential(),dl);

  // Get random number from Landau distribution
  // Float energy_loss = rng.Landau(p.mean,p.sigma);
  Float energy_loss = ThrowLandauHost(p.mpv,4*p.xi,rng.Rndm());

  // Update track
  track.UpdateEnergy(-energy_loss * dl);
}

} // End namespace na63