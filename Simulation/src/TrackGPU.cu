#include "Geometry/BoxCUDA.cuh"
#include "Geometry/SphereCUDA.cuh"
#include "Simulation/TrackGPU.cuh"

namespace na63 {

__device__ inline
bool InsideVolume(const GPUFourVector& position, const VolumePars& volume) {
  const void *specific = (void*)volume.specific;
  if (volume.function_index == SPHERE) {
    return Sphere_Inside(position,specific);
  }
  if (volume.function_index == BOX) {
    return Box_Inside(position,specific);
  }
  return false;
}

__device__
int VolumeQuery(const GPUTrack& track, const VolumePars* volumes,
    const unsigned n_volumes) {

  // Check bounds first
  if (!InsideVolume(track.position,volumes[0])) return -1;

  // Check if still inside current object, if any
  int current = track.volume_index;
  if (current > 0) {
    if (InsideVolume(track.position,volumes[current])) return current;
  }

  // Otherwise check other objects
  for (int i=1;i<n_volumes;i++) {
    if (InsideVolume(track.position,volumes[i])) return i;
  }

  // If none are found, return 0
  return 0;
}

void Boost(GPUTrack& track, const Float bx, const Float by, const Float bz) {

  const Float b2 = pow(bx,2) + pow(by,2) + pow(bz,2);

  const Float gamma = 1.0 / sqrt(1.0 - b2);
  const Float bp = bx*track.position[0] + by*track.position[1] + bz*track.position[2];
  const Float gamma2 = b2 > 0 ? (gamma - 1.0) / b2 : 0.0;

  track.position[0] += gamma2*bp*bx + gamma*bx*track.position[3];
  track.position[1] += gamma2*bp*by + gamma*by*track.position[3];
  track.position[2] += gamma2*bp*bz + gamma*bz*track.position[3];
  track.position[3] = gamma * (track.position[3] + bp);

}


__device__
void Step(GPUTrack& track, const ParticlePars& particle, const Float dl) {
  if (!track.alive) return;
  GPUThreeVector change;
  CUDA_SphericalToCartesian(
    change,
    dl,
    CUDA_CartesianToSpherical_Theta(track.momentum[0],
                                    track.momentum[1],
                                    track.momentum[2]),
    CUDA_CartesianToSpherical_Phi(track.momentum[0],
                                  track.momentum[1])
  );
  track.position[0] += change[0];
  track.position[1] += change[1];
  track.position[2] += change[2];
  track.position[3] += dl * track.momentum[3] / 
      (CUDA_MomentumMagnitude(track.momentum[3],particle.mass) * kC);
}

__device__
void CUDA_UpdateMomentum(GPUFourVector& momentum, const Float change) {
  Float length = CUDA_CartesianToSpherical_R(momentum[0],momentum[1],momentum[2]);
  momentum[0] += change * momentum[0] / length;
  momentum[1] += change * momentum[1] / length;
  momentum[2] += change * momentum[2] / length;
  momentum[3] += change;
}

} // End namespace na63