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
  rng_energy.SetSeed(clock());
}

void Bremsstrahlung(Track& mother, const Material& material,
    const Float dl) {

  if (mother.particle_id != 11) return; // Only treat electrons

  // See if anything happens
  Float chance_to_interact = dl / 10;
  if (rng_interact.Rndm() > chance_to_interact) return;

  // Distribute energy
  Float energy_child = rng_energy.Rndm() * mother.energy();
  mother.UpdateMomentum(-energy_child);

  // Draw angles
  Float theta = rng_theta.Rndm();
  Float phi = rng_phi.Rndm();

  // Spawn child photon
  Float length = CartesianToSpherical_R(mother.momentum[0],mother.momentum[1],mother.momentum[2]);
  FourVector momentum_child(
    energy_child * mother.momentum[0] / length,
    energy_child * mother.momentum[1] / length,
    energy_child * mother.momentum[2] / length,
    energy_child
  );
  Track child(22,mother.position,momentum_child);
  mother.SpawnChild(child);

}

} // End namespace na63