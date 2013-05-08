#ifndef NA63_SIMULATION_BETHEENERGYLOSS_H
#define NA63_SIMULATION_BETHEENERGYLOSS_H

#include <cmath>
#include <iostream>
#include "Geometry/Library.hh"

// K/A = 0.307075 MeV g^-1 cm^2

namespace na63 {

void BetheEnergyLoss(Track& track, const Material& material,
    const Float dl) {

  // Only treat particles with charge
  if (track.charge() == 0) return;

  Float beta_squared = pow(track.beta(),2);
  Float gamma = track.gamma();
  Float gamma_squared = pow(gamma,2);
  Float mass = track.mass();

  // PDG p. 287
  Float T_max = 2 * kElectronMass * beta_squared * gamma_squared /
      (1 + 2 * kElectronMass * gamma / mass + pow(kElectronMass / mass,2));

  // PDG p. 286
  Float mean_energy_loss = material.atomic_number() * pow(track.charge(),2)
      * 0.307075 / beta_squared;
  mean_energy_loss *= 1/2
      * log(2 * kElectronMass * beta_squared * gamma_squared * T_max)
      / material.mean_excitation_potential()
      - beta_squared;

  // Update track
  track.momentum[3] -= mean_energy_loss * dl;
}

} // End namespace na63

#endif /* NA63_SIMULATION_BETHEENERGYLOSS_H */