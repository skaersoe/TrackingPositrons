#include "Geometry/Library.hh"

namespace na63 {

ThreeVector& ThreeVector::operator=(const FourVector& fv) {
  vector[0] = fv[0];
  vector[1] = fv[1];
  vector[2] = fv[2];
  return *this;
}

Float FourVector::length() const {
  return sqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2]);
}

Float Beta(const Float& gamma) {
  if (gamma < 1.0) {
    printf("Received gamma = %g\n",gamma);
    assert(gamma >= 1.0);
  }
  Float beta = sqrt(1-1/pow(gamma,2));
  if (beta == 1) {
    printf("WARNING: float precision exceeded, beta -> 1\n");
    beta = 1 - 1e-6;
  }
  return beta;
}

} // End namespace na63