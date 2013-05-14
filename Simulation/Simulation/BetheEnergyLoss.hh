#ifndef NA63_SIMULATION_BETHEENERGYLOSS_H
#define NA63_SIMULATION_BETHEENERGYLOSS_H

#include "Geometry/Library.hh"
#include "Geometry/Material.hh"
#include "Simulation/Track.hh"

// K/A = 0.307075 MeV g^-1 cm^2

TRandom3 rng;

namespace na63 {

typedef struct {
  Float mean;
  Float sigma;
} LandauParameters;

inline LandauParameters GetSkewedLandauParameters(const Float beta,
    const Float mass, const Float charge, const Float atomic_number,
    const Float mean_excitation_potential, const Float dl) {

  LandauParameters p;

  // Calculate necessary values
  const Float beta_squared = pow(beta,2);
  const Float gamma = Gamma(beta);
  const Float gamma_squared = pow(gamma,2);
  const Float xi = 4.0 * 0.307075 * atomic_number * dl / beta_squared;

  p.mean = xi * (log(2 * kElectronMass * pow(kC,2) * beta_squared
      * gamma_squared / mean_excitation_potential)
      + log(xi/mean_excitation_potential) + 0.200 - beta_squared);
  p.sigma = 4 * xi;

  return p;

}

inline LandauParameters GetBetheLandauParameters(const Float beta,
    const Float mass, const Float charge, const Float atomic_number,
    const Float mean_excitation_potential, const Float dl) {

  LandauParameters p;

  // Grab and calculate necessary values
  const Float beta_squared = pow(beta,2);
  const Float gamma = Gamma(beta);
  const Float gamma_squared = pow(gamma,2);

  // [PDG 27.2.2, eq. 27.4]
  Float T_max = 2.0 * kElectronMass * pow(kC,2) * beta_squared * gamma_squared /
        (1.0 + 2.0 * kElectronMass * gamma / mass + pow(kElectronMass / mass,2));

  // [PDG 27.2.2, eq. 27.3]
  p.mean = pow(charge,2) * atomic_number
      * 0.307075 / beta_squared;
  p.mean *= (0.5 * log((2.0 * kElectronMass * pow(kC,2)
      * beta_squared * gamma_squared * T_max) 
      / pow(mean_excitation_potential,2)) - beta_squared);

  // xi = (K/A) * Z * (x/beta^2)
  // sigma = 4 * xi [PDG 27.2.7]
  p.sigma = 4.0 * 0.307075 * atomic_number * dl / beta_squared;

  return p;
}

void BetheEnergyLoss(Track& track, const Material& material,
    const Float dl) {

  rng.SetSeed((size_t)&track);

  // Only treat particles with charge
  if (track.charge() == 0) return;

  // Get -<dE/dx> and sigma
  LandauParameters p = GetSkewedLandauParameters(
      track.beta(),track.mass(),track.charge(),material.atomic_number(),
      material.mean_excitation_potential(),dl);

  // Get random number from Landau distribution
  Float energy_loss = rng.Landau(p.mean,p.sigma);

  // Update track
  track.momentum[3] -= (energy_loss * dl) / kC;
}

} // End namespace na63

#endif /* NA63_SIMULATION_BETHEENERGYLOSS_H */