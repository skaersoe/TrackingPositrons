#include "Simulation/BremsstrahlungCUDA.cuh"
#include "Simulation/Spawner.cuh"
#include "Geometry/Constants.hh"
#include "Simulation/DeviceGlobalVariables.cuh"
#include "Geometry/LibraryCUDA.cuh"

namespace na63 {

__device__ inline
void CUDA_BremsstrahlungPhoton(GPUTrack& mother,
    const MaterialPars& material, const ParticlePars& particle,
    const Float dl, curandState *rng_state, const int index) {

  // Must have sufficient energy to be considered
  if (mother.momentum[3] <= 2*kElectronMass) return;

  // Must be able to spawn two children
  int child_destination_1 = CanHaveTwoChildren(index);
  if (child_destination_1 == -1) {
    UpdateState(index,WAITING);
    return;
  }
  int child_destination_2 = child_destination_1 - 1;

  // See if anything happens
  Float chance_to_interact = dl / material.radiation_length;
  if (curand_uniform(rng_state) > chance_to_interact) return;

  // Get electron/positron energy and momentum
  Float child_energy = mother.momentum[3] / 2;
  Float child_momentum = sqrt(child_energy*child_energy - kElectronMass*kElectronMass);
  Float phi = kPi * curand_uniform(rng_state);
  Float theta = 2 * kPi * curand_uniform(rng_state);
  GPUThreeVector direction;
  CUDA_SphericalToCartesian(direction,child_momentum,phi,theta);

  // Spawn electron and positron back to back
  GPUTrack *electron = &tracks[child_destination_1];
  electron->particle_id = 11;
  electron->particle_index = electron_index;
  electron->charge = -1;
  electron->momentum[0] = direction[0];
  electron->momentum[1] = direction[1];
  electron->momentum[2] = direction[2];
  electron->momentum[3] = child_energy;
  GPUTrack *positron = &tracks[child_destination_2];
  *positron = *electron;
  positron->charge = 1;
  ThreeVector_Negate(positron->momentum,positron->momentum);

  // Queue new particles for propagation and murder photon
  UpdateState(index,DEAD);
  SpawnChild(mother,child_destination_1);
  SpawnChild(mother,child_destination_2);

}

__device__ inline
void CUDA_BremsstrahlungElectron(GPUTrack& mother,
    const MaterialPars& material, const ParticlePars& particle,
    const Float dl, curandState *rng_state, const int index) {

  // Must be able to spawn a child
  int child_destination = CanHaveChild(index);
  if (child_destination == -1) {
    UpdateState(index,WAITING);
    return;
  }

  // See if anything happens
  Float chance_to_interact = dl / material.radiation_length;
  if (curand_uniform(rng_state) > chance_to_interact) return;

  // Distribute energy between 10% and 30% of electron energy
  Float photon_energy = mother.momentum[3] * (0.10 + 0.20 * curand_uniform(rng_state));

  // Photon direction
  Float mother_phi = CUDA_CartesianToSpherical_Phi(mother.momentum[0],mother.momentum[1]);
  Float mother_theta = CUDA_CartesianToSpherical_Theta(mother.momentum[0],mother.momentum[1],mother.momentum[2]);
  Float mother_gamma = CUDA_Gamma(mother.momentum[3],particle.mass);
  Float phi = mother_phi + mother.charge * asin(1 / mother_gamma);
  GPUThreeVector photon_direction;
  CUDA_SphericalToCartesian(photon_direction,photon_energy,mother_theta,phi);

  // Create new track
  GPUTrack *photon = &tracks[child_destination];
  photon->particle_id = 22;
  photon->particle_index = photon_index;
  ThreeVector_Copy(photon->momentum,photon_direction);
  photon->momentum[3] = photon_energy;

  // Subtract from mother
  FourVector_Subtract(mother.momentum,photon->momentum,mother.momentum);
  
  // Spawn photon
  SpawnChild(mother,child_destination);

}

__device__
void CUDA_Bremsstrahlung(
    GPUTrack& mother,
    const ParticlePars& particle,
    const MaterialPars& material,
    const Float dl,
    curandState *rng_state, 
    const int index) {

  int id = mother.particle_id;
  if (id == 11) CUDA_BremsstrahlungElectron(mother,material,particle,dl,rng_state,index);
  if (id == 22) CUDA_BremsstrahlungPhoton(mother,material,particle,dl,rng_state,index);

}

} // End namespace na63