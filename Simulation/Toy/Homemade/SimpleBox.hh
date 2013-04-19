#ifndef SIMPLEBOX_H
#define SIMPLEBOX_H

#include "Track.hh"
#include "Types.hh"

class SimpleBox {

public:
  SimpleBox(ThreeVector c, ThreeVector s)
      : center(c), sides(s) {}
  bool Inside(ThreeVector *p) {
    return p->x <= center.x + sides.x/2 &&
           p->x >= center.x - sides.x/2 &&
           p->y <= center.y + sides.y/2 &&
           p->y >= center.y - sides.y/2 &&
           p->z <= center.z + sides.z/2 &&
           p->z >= center.z - sides.z/2;
  }
  bool Inside(Track *t) {
    return Inside(&t->position);
  }

private:
  ThreeVector center;
  ThreeVector sides;

};

#endif