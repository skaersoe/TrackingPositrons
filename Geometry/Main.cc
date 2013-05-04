#include <iostream>
#include "Geometry/Library.hh"
#include "Geometry/Geometry.hh"
#include "Geometry/Box.hh"
#include "Geometry/Sphere.hh"
#include "Simulation/Simulator.hh"

using namespace na63;

int main(void) {

  // Create some geometry
  Geometry geometry;
  geometry.AddMaterial(Material("iron",1,1));
  Box box = Box("iron",ThreeVector(0,0,0),ThreeVector(2e2,2e2,2e2));
  Sphere sphere = Sphere("iron",ThreeVector(0,0,0),2e2);
  geometry.AddVolume(box);
  geometry.AddVolume(sphere);

  Track t = Track(11,FourVector(1,1,1,1),FourVector(0,0,0,0));

  std::cout << "Track inside: " << box.Inside(t.position) << std::endl;
  std::cout << "Track inside: " << sphere.Inside(t.position) << std::endl;

  return 0;
}