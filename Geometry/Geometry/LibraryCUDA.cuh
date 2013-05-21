#ifndef NA63_LIBRARYCUDA_CUH
#define NA63_LIBRARYCUDA_CUH

#include "Geometry/Library.hh"
#include <math.h>

namespace na63 {

#define MAX_INT_VALUE 2147483647

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

__device__ inline
Float CUDA_CartesianToSpherical_R(const Float x, const Float y, const Float z) {
  return sqrt(pow(x,2) + pow(y,2) + pow(z,2));
}

__device__ inline
Float CUDA_CartesianToSpherical_Theta(const Float x, const Float y, const Float z) {
  return acos(z/CUDA_CartesianToSpherical_R(x,y,z));
}

__device__ inline
Float CUDA_CartesianToSpherical_Phi(const Float x, const Float y) {
  return atan(y/x);
}

__device__ inline
Float CUDA_Gamma(const Float& beta) {
  return 1/sqrt(1-pow(beta,2));
}

__device__ inline
Float CUDA_Gamma(const Float& energy, const Float& mass) {
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
Float CUDA_Beta(const Float& energy, const Float& mass) {
  Float energy_squared = pow(energy,2);
  return sqrt((energy_squared - pow(mass,2)) / energy_squared);
}

} // End namespace na63

#endif /* NA63_LIBRARYCUDA_CUH */