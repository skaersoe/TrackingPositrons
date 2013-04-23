#include <iostream>
#include <time.h>
#include <cassert>

#include "Simulation/Simulator.hh"
#include "Geometry/Geometry.hh"

using namespace na63;

int main(int argc,char *argv[]) {

  // Set some default values
  SimulatorDevice device = GPU;
  bool debug = false;
  unsigned N = 0;

  // Create geometry
  Geometry geometry;
  geometry.AddMaterial(Material("vacuum",0.0));
  geometry.AddMaterial(Material("solid",1.0));
  geometry.AddVolume(Sphere("solid",{0,0,0},5));
  geometry.AddVolume(Sphere("vacuum",{0,0,0},100));
  geometry.AddParticle(Particle("electron",11,-1,0.510998910));

  // Create Simulator object
  Simulator *sim = new Simulator(&geometry);

  // Parse input
  if (argc < 2) {
    std::cerr << "Missing input: number of particles to simulate." << std::endl;
    return -1;
  }
  sscanf(argv[1],"%d",&N);
  for (int i=2;i<argc;i++) {
    std::string token(argv[i]);
    if (token == "-CPU")
      device = CPU;
    else if (token == "-debug")
      debug = true;
    else {
      std::cout << "Unrecognized argument: " << token << std::endl;
       return -1;
    }
  } // Finished parsing input arguments

  sim->pars.device = device;
  sim->pars.debug = debug;
  sim->pars.N = N;

  std::cout << "Simulator starting." << std::endl;

  // Populate
  clock_t timer = clock();
  sim->GenerateTracks();
  timer = clock() - timer;
  std::cout << "Populated in " << (float)timer/CLOCKS_PER_SEC << " seconds." << std::endl;

  // Propagate
  timer = clock();
  sim->Propagate();
  timer = clock() - timer;
  std::cout << "Propagated in " << (float)timer/CLOCKS_PER_SEC << " seconds." << std::endl;

  delete sim;

  std::cout << "Simulator exiting." << std::endl;
  return 0;
}