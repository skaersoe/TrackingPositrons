#include "Simulation/AngularDistributionCUDA.cuh"
#include "Geometry/LibraryCUDA.cuh"
#include "Geometry/Constants.hh"

#include <curand_kernel.h>

namespace na63 {

__device__
void CUDA_ModifiedTsai_SampleDirection(GPUThreeVector& local_direction,
    const GPUTrack* track, const Float mass, curandState *rng_state) {

  // Sample gamma angle (Z - axis along the parent particle).
  // Universal distribution suggested by L. Urban (Geant3 manual (1993) 
  // Phys211) derived from Tsai distribution (Rev Mod Phys 49,421(1977))
  
  Float u_max = 2.0 * (1.0 + (track->momentum[3] - mass)/kElectronMass);   

  Float a1     = 0.625;
  Float a2     = 1.875;
  Float border = 0.25;
  Float u;

  do {
    u = - log(curand_uniform(rng_state)*curand_uniform(rng_state));
    if (border > curand_uniform(rng_state)) {
      u /= a1;
    } else {
      u /= a2;
    }

  } while (u > u_max);

  Float cost = 1.0 - 2*u*u/(u_max*u_max);
  Float sint = sqrt((1 - cost)*(1 + cost));
  Float phi  = 2*kPi*curand_uniform(rng_state); 

  ThreeVector_Set(local_direction,sint*cos(phi), sint*sin(phi), cost);
  GPUThreeVector norm;
  ThreeVector_Normalized(norm,track->momentum);
  ThreeVector_Rotate(local_direction,norm);

}

} // End namespace na63