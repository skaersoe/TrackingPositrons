#include "Simulation/Bremsstrahlung.hh"
#include "Geometry/Library.hh"

namespace na63 {

TRandom3 rng_interact;
TRandom3 rng_energy;
TRandom3 rng_theta;
TRandom3 rng_phi;

void InitializeBremsstrahlung() {
  rng_interact.SetSeed(clock());
  rng_energy.SetSeed(clock());
  rng_theta.SetSeed(clock());
  rng_phi.SetSeed(clock());
}

inline void BremsstrahlungPhoton(Track& mother, const Material& material,
    const Float dl) {

  // Must have enough energy to be considered
  if (mother.energy() <= 2*kElectronMass) return;

  // See if anything happens
  Float chance_to_interact = dl / material.radiation_length();
  if (rng_interact.Rndm() > chance_to_interact) return;

  // Get electron/positron energy and momentum
  Float child_energy = mother.energy() / 2;
  Float child_momentum = sqrt(child_energy*child_energy - kElectronMass*kElectronMass);
  Float phi = kPi * rng_phi.Rndm();
  Float theta = 2 * kPi * rng_theta.Rndm();
  ThreeVector direction = SphericalToCartesian(child_momentum,phi,theta);

  // Spawn electron and positron back to back and boost them
  Track electron(11,-1,mother.position,FourVector( direction[0], direction[1], direction[2],child_energy));
  Track positron(11,-1,mother.position,FourVector(-direction[0],-direction[1],-direction[2],child_energy));
  ThreeVector beta = mother.beta_vector();
  electron.Boost(beta);
  positron.Boost(beta);

  // Queue new particles for propagation and murder photon
  mother.SpawnChild(electron);
  mother.SpawnChild(positron);
  mother.Kill();

}

inline void BremsstrahlungElectron(Track& mother, const Material& material,
    const Float dl) {

  // See if anything happens
  Float chance_to_interact = dl / material.radiation_length();
  if (rng_interact.Rndm() > chance_to_interact) return;

  // Distribute energy from 10% to 30% of electron energy
  Float photon_energy = mother.energy() * (0.1 + 0.2 * rng_energy.Rndm());

  // Photon direction
  Float phi = mother.phi() + mother.charge() * asin(1 / mother.gamma());
  ThreeVector photon_direction = SphericalToCartesian(photon_energy,mother.theta(),phi);
  Track photon(22,0,mother.position,FourVector(photon_direction[0],photon_direction[1],photon_direction[2],photon_energy));

  // Subtract from mother
  mother.momentum -= photon.momentum;
  
  // Spawn photon
  mother.SpawnChild(photon);

}

void Bremsstrahlung(Track& mother, const Material& material,
    const Float dl) {

  if (mother.particle_id == 11) BremsstrahlungElectron(mother,material,dl);
  if (mother.particle_id == 22) BremsstrahlungPhoton(mother,material,dl);

}

} // End namespace na63