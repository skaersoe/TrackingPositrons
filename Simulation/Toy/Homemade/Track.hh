#ifndef TRACK_H
#define TRACK_H

#include <cmath>
#include <iostream>
#include "Types.hh"

class SimpleBox;

class Track {

public:
  Track(ThreeVector pos, FourMomentum mom, float cha)
      : position(pos), momentum(mom) {
    charge = cha;
  }
  float x() { return position.x; }
  float y() { return position.y; }
  float z() { return position.z; }
  float px() { return momentum.px(); }
  float py() { return momentum.py(); }
  float pz() { return momentum.pz(); }
  FourMomentum P() { return momentum; }
  float q() { return charge; }
  float E() { return momentum.E(); }
  void timestep(float dt) {
    position.x += dt * px();
    position.y += dt * py();
    position.z += dt * pz();
  }

private:
  ThreeVector position;
  FourMomentum momentum;
  float charge;
  friend class SimpleBox;

};

#endif