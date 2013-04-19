#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "SimpleBox.hh"
#include <vector>

class Geometry {
public:
  Geometry(SimpleBox b) : bounds(b) {}
  void AddBox(SimpleBox b) {
    boxes.push_back(b);
  }
  bool InBounds(Track *track) {
    return bounds.Inside(track);
  }
  SimpleBox* BoundaryDetection(Track *track) {
    for (int i=0;i<boxes.size();i++) {
      if (boxes[i].Inside(track)) return &boxes[i];
    }
    return NULL;
  }

private:
  SimpleBox bounds;
  std::vector<SimpleBox> boxes;

};

#endif