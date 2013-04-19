#ifndef TYPES_H
#define TYPES_H

class ThreeVector {
public:
  ThreeVector(float xi, float yi, float zi) {
    x = xi;
    y = yi;
    z = zi;
  }
  ThreeVector(const ThreeVector& v) {
    x = v.x;
    y = v.y;
    z = v.z;
  }
  float x, y, z;
  ThreeVector& operator= (const ThreeVector& rhs) {
    x = rhs.x;
    y = rhs.y;
    z = rhs.z;
    return *this;
  }
};

class FourMomentum {
public:
  FourMomentum(ThreeVector m, float e) : momentum(m) {
    energy = e;
  }
  float px() { return momentum.x; }
  float py() { return momentum.y; }
  float pz() { return momentum.z; }
  ThreeVector p() { return momentum; };
  float E() { return energy; }
  FourMomentum& operator= (const FourMomentum& rhs) {
    momentum = rhs.momentum;
    energy = rhs.energy;
    return *this;
  }
private:
  ThreeVector momentum;
  float energy;
};

#endif