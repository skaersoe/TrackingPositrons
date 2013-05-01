#ifndef TRACK_H
#define TRACK_H

#include <cmath>
#include <iostream>
#include "Library.hh"

class SimpleBox;

class Track {

public:
  FourVector position;
  FourVector momentum;
  Float charge;
  Float mass;

  Track(FourVector pos, FourVector mom, Float cha, Float ma)
      : position(pos), momentum(mom) {
    charge = cha;
    mass = ma;
  }

  inline Float x() const { return position[0]; }
  inline Float y() const { return position[1]; }
  inline Float z() const { return position[2]; }
  inline Float t() const { return position[3]; }
  inline Float px() const { return momentum[0]; }
  inline Float py() const { return momentum[1]; }
  inline Float pz() const { return momentum[2]; }
  inline Float E() const { return momentum[3]; }
  inline Float q() const { return charge; }
  inline Float m() const { return mass; }
  inline Float p() const {
    return sqrt(pow(E(),2) - pow(m(),2));
  }
  inline Float beta() const {
    Float E_squared = pow(E(),2);
    return sqrt((E_squared-pow(m(),2))/E_squared);
  }
  inline Float gamma() const {
    return E() / m();
  }
  ostream& operator<< (ostream& os) {
    os << "(" << position[0]
       << "," << position[1]
       << "," << position[2]
       << "," << position[3] << "), "
       << "(" << momentum[0]
       << "," << momentum[1]
       << "," << momentum[2]
       << "," << momentum[3] << ")";
    return os;
  }

  void Step(Float dl);

};

#endif