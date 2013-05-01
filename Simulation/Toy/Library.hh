#ifndef TYPES_H
#define TYPES_H

#include "TLorentzVector.h"

typedef double Float;
typedef TLorentzVector FourVector;

class ThreeVector {
private:
  Float vector[3];
public:
  Float& operator[] (const int i) {
    return vector[i];
  }
  friend ostream& operator<< (ostream& os, const ThreeVector& tv) {
    os << "(" << tv.vector[0]
       << "," << tv.vector[1]
       << "," << tv.vector[2] << ")";
    return os;
  }
  ThreeVector(Float a, Float b, Float c) {
    vector[0] = a;
    vector[1] = b;
    vector[2] = c;
  }
};

const Float c = 299792458;
const Float pi = 3.14159265;

inline ThreeVector SphericalToCartesian(Float r, Float theta, Float phi) {
  return ThreeVector(r * sin(theta) * cos(phi), // x
                     r * sin(theta) * sin(phi), // y
                     r * cos(theta));           // z
}

inline Float Coulomb(Float e) {
  return 1.602176565e-19 * e;
}

inline Float JouleToGev(Float J) {
  return 6.24150934e9 * J;
}

#endif