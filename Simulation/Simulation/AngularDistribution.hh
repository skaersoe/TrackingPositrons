#ifndef NA63_SIMULATION_ANGULARDISTRIBUTION_H
#define NA63_SIMULATION_ANGULARDISTRIBUTION_H

#include <cmath>

#include "Geometry/Library.hh"
#include "Geometry/Constants.hh"

namespace na63 {

/** Generates random angular distributions for bremsstrahlung */

ThreeVector ModifiedTsai_SampleDirection(const Track& track) {

  // Sample gamma angle (Z - axis along the parent particle).
  // Universal distribution suggested by L. Urban (Geant3 manual (1993) 
  // Phys211) derived from Tsai distribution (Rev Mod Phys 49,421(1977))
  
  ThreeVector local_direction;
  Float u_max = 2.0 * (1.0 + track.kinetic_energy()/kElectronMass);   

  const Float a1     = 0.625;
  const Float a2     = 1.875;
  const Float border = 0.25;
  Float u;

  do {
    u = - log(gRandom->Uniform()*gRandom->Uniform());
    if (border > gRandom->Uniform()) {
      u /= a1;
    } else {
      u /= a2;
    }

  } while (u > u_max);

  Float cost = 1.0 - 2*u*u/(u_max*u_max);
  Float sint = sqrt((1 - cost)*(1 + cost));
  Float phi  = 2*kPi*gRandom->Uniform(); 

  local_direction.Set(sint*cos(phi), sint*sin(phi), cost);
  local_direction.Rotate(track.momentum.Normalized());

  return local_direction;

}

} // End namespace na63

#endif /* NA63_SIMULATION_ANGULARDISTRIBUTION_H */