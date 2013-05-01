#include <iostream>
#include <time.h>
#include <cassert>

#include "Simulation/Simulator.hh"
#include "Geometry/Geometry.hh"

using namespace na63;

int main(int argc,char *argv[]) {

  // Create geometry
  Geometry geometry;
  geometry.AddMaterial(Material("iron",26,286));
  geometry.AddParticle(Particle("electron",11,-1,0.510998910));

  // Create Simulator object
  Simulator sim = Simulator(&geometry);

  // Set some default values
  sim.device = GPU;
  sim.debug = false;
  unsigned N = 0;

  // Parse input
  if (argc < 2) {
    std::cerr << "Missing input: number of particles to simulate." << std::endl;
    return -1;
  }
  sscanf(argv[1],"%d",&N);
  for (int i=2;i<argc;i++) {
    std::string token(argv[i]);
    if (token == "-CPU")
      sim.device = CPU;
    else if (token == "-debug")
      sim.debug = true;
    else {
      std::cout << "Unrecognized argument: " << token << std::endl;
      return -1;
    }
  } // Finished parsing input arguments

  std::cout << "Simulator starting." << std::endl;

  // Populate
  Track t[5] = Track(11,FourVector(1,2,3,4),FourVector(5,6,7,8));
  sim.SetTracks(&t[0],5);

  // Propagate
  clock_t timer = clock();
  sim.Propagate();
  timer = clock() - timer;
  std::cout << "Propagated in " << (float)timer/CLOCKS_PER_SEC << " seconds." << std::endl;

  std::cout << "Simulator exiting." << std::endl;
  return 0;
}