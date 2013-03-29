#include <iostream>
#include <time.h>

#include "Simulation/Simulator.hh"
#include "Simulation/Track.hh"

using namespace na63;

int main(int argc,char *argv[]) {

  // Set some default values
  SimulatorDevice device = GPU;
  bool debug = false;
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
      device = CPU;
    else if (token == "-debug")
      debug = true;
    else {
      std::cout << "Unrecognized argument: " << token << std::endl;
       return -1;
    }
  } // Finished parsing input arguments

  std::cout << "Simulator starting." << std::endl;

  // Create Simulator object
  Simulator *sim = new Simulator(device, debug, N);

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