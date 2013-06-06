#include <iostream>
#include "Geometry/SimpleBox.hh"

using namespace na63;

int main(void) {
  FourVector pos(0,-1,0,0);
  SimpleBox b("derp",ThreeVector(2e3,0,0),ThreeVector(4e3,4e3,4e3));
  std::cout << b.Inside(pos) << std::endl;
  return 0;
}