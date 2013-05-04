#ifndef NA63_LIBRARYCUDA_CUH
#define NA63_LIBRARYCUDA_CUH

#include "Geometry/Library.hh"
#include <math.h>

namespace na63 {

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

} // End namespace na63

#endif /* NA63_LIBRARYCUDA_CUH */