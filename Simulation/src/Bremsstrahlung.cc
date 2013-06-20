#include "Simulation/Bremsstrahlung.hh"
#include "Simulation/AngularDistribution.hh"
#include "Geometry/Material.hh"

namespace na63 {

inline void BremsstrahlungPhoton(Track& mother, const Material& material,
    const Float dl) {

  // Must have enough energy to be considered
  if (mother.energy() <= 2*kElectronMass) {
    //mother.Kill();
    return;
  }

  // See if anything happens
  Float chance_to_interact = dl / material.radiation_length();
  if (gRandom->Uniform() > chance_to_interact) return;

  // Get electron/positron energy and momentum
  Float child_energy = mother.energy() / 2;
  Float child_momentum = sqrt(child_energy*child_energy - kElectronMass*kElectronMass);
  Float phi = kPi * gRandom->Uniform();
  Float theta = 2 * kPi * gRandom->Uniform();
  ThreeVector direction = SphericalToCartesian(child_momentum,phi,theta);
  direction.Extend(sqrt(child_momentum));

  // Spawn electron and positron back to back and boost them
  Track electron(11,-1,mother.position,FourVector( direction[0], direction[1], direction[2],child_energy));
  Track positron(11, 1,mother.position,FourVector(-direction[0],-direction[1],-direction[2],child_energy));

  // Queue new particles for propagation and murder photon
  mother.SpawnChild(electron);
  mother.SpawnChild(positron);
  mother.Kill();

}

inline void BremsstrahlungElectron(Track& mother, const Material& material,
    const Float dl) {

  // See if anything happens
  Float chance_to_interact = dl / material.radiation_length();
  if (gRandom->Uniform() > chance_to_interact) return;

  // Distribute energy from 10% to 30% of electron energy
  Float photon_energy = mother.energy() * (0.1 + 0.2 * gRandom->Uniform());

  // Photon momentum
  ThreeVector photon_momentum = ModifiedTsai_SampleDirection(mother);
  photon_momentum.Extend(sqrt(photon_energy));
  Track photon(22,0,mother.position,FourVector(photon_momentum[0],photon_momentum[1],photon_momentum[2],photon_energy));

  // Subtract from mother
  mother.UpdateMomentum(photon.momentum);
  
  // Spawn photon
  mother.SpawnChild(photon);

}

void MyBremsstrahlung(Track& mother, const Material& material,
    const Float dl) {

  if (mother.particle_id == 11) BremsstrahlungElectron(mother,material,dl);
  if (mother.particle_id == 22) BremsstrahlungPhoton(mother,material,dl);

}

} // End namespace na63