#include <iostream>
#include "Geometry/Library.hh"
#include "Geometry/Geometry.hh"
#include "Geometry/Sphere.hh"
#include "Simulation/Simulator.hh"

using namespace na63;

int main(void) {

  // Create some geometry
  Geometry geometry;
  geometry.AddMaterial(Material("iron",0,0));
  geometry.AddVolume(Sphere("solid",{0,0,0},5));
  geometry.AddVolume(Sphere("vacuum",{0,0,0},100));
  geometry.AddParticle(Particle("electron",11,-1,0.510998910));
  geometry.PrintContent();



  return 0;
}