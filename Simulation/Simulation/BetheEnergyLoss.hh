#ifndef NA63_SIMULATION_BETHEENERGYLOSS_H
#define NA63_SIMULATION_BETHEENERGYLOSS_H

#include "Geometry/Library.hh"
#include "Geometry/Constants.hh"
#include "Geometry/Material.hh"
#include "Simulation/Track.hh"

// K/A = 0.307075 MeV g^-1 cm^2

namespace na63 {

typedef struct {
  Float mpv;
  Float xi;
} LandauParameters;

inline LandauParameters LandauEnergyLossParameters(const Float beta,
    const Float mass, const Float atomic_number,
    const Float mean_excitation_potential, const Float dl) {

  LandauParameters p;

  // Calculate necessary values
  Float beta_squared = pow(beta,2);
  Float gamma = Gamma(beta);
  Float gamma_squared = pow(gamma,2);
  p.xi = 0.5 * 0.307075 * atomic_number * dl / beta_squared;

  p.mpv = p.xi * (log(2 * mass * beta_squared
      * gamma_squared / mean_excitation_potential)
      + log(p.xi/mean_excitation_potential) + 0.200 - beta_squared);

  return p;

}

// inline LandauParameters GetBetheLandauParameters(const Float beta,
//     const Float mass, const Float charge, const Float atomic_number,
//     const Float mean_excitation_potential, const Float dl) {

//   LandauParameters p;

//   // Grab and calculate necessary values
//   const Float beta_squared = pow(beta,2);
//   const Float gamma = Gamma(beta);
//   const Float gamma_squared = pow(gamma,2);

//   // [PDG 27.2.2, eq. 27.4]
//   Float T_max = 2.0 * kElectronMass * beta_squared * gamma_squared /
//         (1.0 + 2.0 * kElectronMass * gamma / mass + pow(kElectronMass / mass,2));

//   // [PDG 27.2.2, eq. 27.3]
//   p.mean = pow(charge,2) * atomic_number
//       * 0.307075 / beta_squared;
//   p.mean *= (0.5 * log((2.0 * kElectronMass
//       * beta_squared * gamma_squared * T_max) 
//       / pow(mean_excitation_potential,2)) - beta_squared);
//   printf("%f, %f\n",pow(charge,2) * atomic_number
//       * 0.307075 / beta_squared,(0.5 * log((2.0 * kElectronMass/* * pow(kC,2)*/
//       * beta_squared * gamma_squared * T_max) 
//       / pow(mean_excitation_potential,2)) - beta_squared));

//   // xi = (K/A) * Z * (x/beta^2)
//   // sigma = 4 * xi [PDG 27.2.7]
//   p.sigma = 2.0 * 0.307075 * atomic_number * dl / beta_squared;

//   return p;
// }

inline Float Bethe_dEdx(const Float beta,
    const Float mass, const Float charge, const Float atomic_number,
    const Float mean_excitation_potential, const Float density,
    const Float atomic_weight, const Float dl) {

  Float beta_squared = beta*beta;
  Float betagamma_squared = pow(Gamma(beta),2) * beta_squared;

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

inline Float BetheElectron_dEdx(const Float density,
    const Float atomic_number, const Float beta, const Float kinetic_energy,
    const Float mean_excitation_potential, const Float charge) {

  // Kinetic energy in terms of electron masses
  Float tau = kinetic_energy / kElectronMass;
  Float beta_squared = beta*beta;

  Float F_tau_electron = 1.0 - beta_squared
      + (tau*tau / 8.0 - log(2.0) * (2.0*kElectronRadius + 1.0)) / pow(tau + 1.0,2);

  Float F_tau_positron = 2 * log(2) - beta_squared / 12
      * (23 + 14/(tau + 2) + 10/pow(tau + 2,2) + 4/pow(tau+2,3));

  Float mean = 2.0 * kPi * kAvogadro * kElectronRadiusSquared * kElectronMass
      * density * atomic_number / (beta*beta);
  mean *= log((tau*tau * (tau + 2.0))
      / (2.0 * pow(1e-6 * mean_excitation_potential / kElectronMass,2)))
      + ((charge < 0) ? F_tau_electron : F_tau_positron);

  return mean;

}

void BetheEnergyLoss(Track& track, const Material& material,
    const Float dl);

} // End namespace na63

#endif /* NA63_SIMULATION_BETHEENERGYLOSS_H */