#ifndef SIMPLEBOX_H
#define SIMPLEBOX_H

#include <iostream>
#include "Track.hh"
#include "Library.hh"

class SimpleBox {

public:
  SimpleBox(ThreeVector c, ThreeVector s)
      : center(c), sides(s) {}
  bool Inside(Float x, Float y, Float z) {
    return x <= center[0] + sides[0]/2 &&
           x >= center[0] - sides[0]/2 &&
           y <= center[1] + sides[1]/2 &&
           y >= center[1] - sides[1]/2 &&
           z <= center[2] + sides[2]/2 &&
           z >= center[2] - sides[2]/2;
  }
  bool Inside(Track *t) {
    return Inside(t->position[0],
                  t->position[1],
                  t->position[2]);
  }

private:
  ThreeVector center;
  ThreeVector sides;

};

#endif