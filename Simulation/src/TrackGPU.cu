#include "Simulation/TrackGPU.cuh"

namespace na63 {

__device__
int VolumeQuery(const GPUFourVector& position, const VolumePars* volumes,
    const InsideFunction* functions, const int volume_index) {
  
  //const VolumePars *v = &volumes[0];
  //InsideFunction f = functions[v->function_index];
  //if (!f(position,(void*)v->specific)) return -1;

  return 0;
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

} // End namespace na63