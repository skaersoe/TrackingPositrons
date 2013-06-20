#ifndef NA63_LIBRARYCUDA_CUH
#define NA63_LIBRARYCUDA_CUH

#include "Geometry/Library.hh"
#include "Simulation/PropagateGPU.cuh"
#include "Simulation/DeviceGlobalVariables.cuh"
#include <math.h>

namespace na63 {

#define MAX_INT_VALUE 2147483647

template <typename A, typename B>
__device__ inline
void ThreeVector_Copy(
    A& a, const B& b) {
  a[0] = b[0];
  a[1] = b[1];
  a[2] = b[2];
}

template <typename A, typename B>
__device__ inline
void FourVector_Copy(
    A& a, const B& b) {
  a[0] = b[0];
  a[1] = b[1];
  a[2] = b[2];
  a[3] = b[3];
}

template <typename A, typename B, typename C>
__device__ inline
void ThreeVector_Add(
    const A& a, const B& b, C& c) {
  c[0] = a[0] + b[0];
  c[1] = a[1] + b[1];
  c[2] = a[2] + b[2];
}

template <typename A, typename B, typename C>
__device__ inline
void ThreeVector_Subtract(
    const A& a, const B& b, C& c) {
  c[0] = a[0] - b[0];
  c[1] = a[1] - b[1];
  c[2] = a[2] - b[2];
}

template <typename A, typename B, typename C>
__device__ inline
void FourVector_Add(
    const A& a, const B& b, C& c) {
  c[0] = a[0] + b[0];
  c[1] = a[1] + b[1];
  c[2] = a[2] + b[2];
  c[3] = a[3] + b[3];
}

template <typename A, typename B, typename C>
__device__ inline
void FourVector_Subtract(
    const A& a, const B& b, C& c) {
  c[0] = a[0] - b[0];
  c[1] = a[1] - b[1];
  c[2] = a[2] - b[2];
  c[3] = a[3] - b[3];
}

template <typename A, typename B>
__device__ inline
Float ThreeVector_Dot(
    const A& a, const B& b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

template <typename A, typename B>
__device__ inline
void ThreeVector_Negate(const A& a, B& b) {
  b[0] = -a[0];
  b[1] = -a[1];
  b[2] = -a[2];
}

template <typename A>
__device__ inline
Float ThreeVector_Length(const A& a) {
  // printf("sqrt(%f^2 + %f^2 + %f^2) = %f\n",a[0],a[1],a[2],sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]));
  return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}

template <typename A>
__device__ inline
void ThreeVector_Extend(A& a, const Float l) {
  a[0] *= l;
  a[1] *= l;
  a[2] *= l;
}

template <typename A>
__device__ inline
void ThreeVector_Normalize(A& a) {
  Float l = 1.0/ThreeVector_Length(a);
  ThreeVector_Extend(a,l);
}

template <typename A, typename B>
__device__ inline
void ThreeVector_Normalized(A& a, const B& b) {
  Float l = 1.0/ThreeVector_Length(b);
  ThreeVector_Copy(a,b);
  ThreeVector_Extend(a,l);
}

__device__ inline
void ThreeVector_Rotate(GPUThreeVector& tgt,
    const GPUThreeVector& unit_vector) {

  Float u1 = unit_vector[0];
  Float u2 = unit_vector[1];
  Float u3 = unit_vector[2];
  Float up = u1*u1 + u2*u2;

  if (up) {
    up = sqrt(up);
    Float px = tgt[0], py = tgt[1], pz = tgt[2];
    tgt[0] = (u1*u3*px - u2*py + u1*up*pz)/up;
    tgt[1] = (u2*u3*px + u1*py + u2*up*pz)/up;
    tgt[2] = (u3*u3*px -    px + u3*up*pz)/up;
  } else {
    if (u3 < 0.0) {
      tgt[0] = -tgt[0]; tgt[2] = -tgt[2];
    }
  }

}

__device__ inline
void ThreeVector_Set(GPUThreeVector& tgt, const Float a, const Float b,
    const Float c) {
  tgt[0] = a;
  tgt[1] = b;
  tgt[2] = c;
}

__device__ inline
void FourVector_Set(GPUFourVector& tgt, const GPUThreeVector& a, const Float b) {
  tgt[0] = a[0];
  tgt[1] = a[1];
  tgt[2] = a[2];
  tgt[3] = b;
}

__device__ inline
Float CUDA_CartesianToSpherical_R(const Float x, const Float y, const Float z) {
  return sqrt(pow(x,2) + pow(y,2) + pow(z,2));
}

__device__ inline
Float CUDA_CartesianToSpherical_Theta(const Float x, const Float y, const Float z) {
  return (!x && !y && !z) ? 0.0 : acos(z/CUDA_CartesianToSpherical_R(x,y,z));
}

__device__ inline
Float CUDA_CartesianToSpherical_Phi(const Float x, const Float y) {
  return (!x && !y) ? 0.0 : atan(y/x);
}

__device__ inline
Float CUDA_Gamma(const Float& beta) {
  // Float gamma = 1/sqrt(1-beta*beta);
  // if (gamma < 1.0) {
  //   printf("Gamma less than one because of beta = %g\n",beta);
  // }
  return 1/sqrt(1-beta*beta);
}

__device__ inline
Float CUDA_Gamma(const Float& energy, const Float& mass) {
  // Float gamma = energy / mass;
  // if (gamma != gamma) {
  //   printf("CUDA_Gamma calculated NaN from %g / %g\n",energy,mass);
  // }
  // if (gamma < 1.0) {
  //   printf("CUDA_Gamma is less than one from %g / %g\n",energy,mass);
  // }
  return energy / mass;
}

__device__ inline
void CUDA_SphericalToCartesian(GPUThreeVector& tv,
    const Float r, const Float theta, const Float phi) {
  tv[0] = r * sin(theta) * cos(phi);
  tv[1] = r * sin(theta) * sin(phi);
  tv[2] = r * cos(theta);
}

__device__ inline
Float CUDA_MomentumMagnitude(const Float& energy, const Float& mass) {
  return sqrt(pow(energy,2) - pow(mass,2));
} 

__device__ inline
Float CUDA_Beta(const Float& gamma) {
  Float beta = sqrt(1-1/(gamma*gamma));
  if (beta == 1) {
    // printf("WARNING: float precision exceeded, beta -> 1\n");
    beta = 1 - 1e-6;
  }
  // if (beta != beta) {
  //   printf("CUDA_Beta calculated NaN from sqrt(1-1/(%g^2)\n",gamma);
  // }
  return beta;
}

__device__ inline
void UpdateState(const int index, const TrackState state) {
  tracks[index].state = state;
  int key = 0;
  if (state == DEAD) {
    key = TRACK_KEY_DEAD;
  }
  else if (state == WAITING) {
    key = TRACK_KEY_WAITING;
  }
  else if (state == FREE) {
    key = TRACK_KEY_FREE;
  }
  keys[index] = key;
}

} // End namespace na63

#endif /* NA63_LIBRARYCUDA_CUH */