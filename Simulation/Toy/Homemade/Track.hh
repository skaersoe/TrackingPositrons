#ifndef TRACK_H
#define TRACK_H

#include <cmath>
#include <iostream>
#include "TLorentzVector.h"

class SimpleBox;

class Track {

public:
  TLorentzVector position;
  TLorentzVector momentum;
  float charge;

  Track(TLorentzVector pos, TLorentzVector mom, float cha)
      : position(pos), momentum(mom) {
    charge = cha;
  }

  double x() { return position[0]; }
  double y() { return position[1]; }
  double z() { return position[2]; }
  double t() { return position[3]; }
  double px() { return momentum[0]; }
  double py() { return momentum[1]; }
  double pz() { return momentum[2]; }
  double E() { return momentum[3]; }
  double q() { return charge; }

};

#endif