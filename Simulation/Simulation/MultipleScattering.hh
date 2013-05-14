#ifndef NA63_SIMULATION_MULTIPLESCATTERING_H
#define NA63_SIMULATION_MULTIPLESCATTERING_H

#include "Geometry/Library.hh"
#include "Geometry/Material.hh"
#include "Simulation/Track.hh"

TRandom3 rng1;
TRandom3 rng2;

namespace na63 {
  
void MultipleScattering(Track& track, const Material& material,
    const Float dl) {

  // Get two independent random numbers
  rng1.SetSeed((size_t)&track);
  rng2.SetSeed((size_t)&track+42);
  Float z1 = rng1.Gaus(0,1);
  Float z2 = rng2.Gaus(0,1);

  // [PDG 27.4.1, eq. 27.24]
  Float atomic_number = material.atomic_number();
  Float X_0 = 716.4 / (atomic_number * (atomic_number + 1.0)
      * log(287.0/sqrt(atomic_number)));

  // [PDG 27.3, eq. 27.14]
  // theta_0 = 13.6 MeV / (beta * c * p) * z * sqrt(x/X_0)
  //           * [1 + 0.038 * ln(x/X_0)]
  Float theta_0 = 13.6 / (track.beta() * kC * track.momentum_magnitude())
      * track.charge() * sqrt(dl / X_0) * (1.0 + 0.038 * log(dl / X_0));

  // [PDG 27.3, eq. 27.20]
  Float y_plane = z1 * dl * theta_0 * 0.2887 + 0.5 * z2 * dl * theta_0;
  // [PDG 27.3, eq. 27.21]
  Float theta_plane = z2 * theta_0;

}

} // End namespace na63

#endif /* NA63_SIMULATION_MULTIPLESCATTERING_H */