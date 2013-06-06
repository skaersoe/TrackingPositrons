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

} // End namespace na63