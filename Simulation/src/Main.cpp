#include <iostream>
#include <time.h>

#include "Simulation/Simulator.h"
#include "Simulation/Particle.h"
#include "Simulation/Arguments.h"

using namespace NA63;

int main(int argc,char *argv[]) {

  /* Set some default values. N must be specified by user */
  simulator_args_t args = { .device = GPU, .debug = false, .render = false, .N = 0 };
  /* Parse input */
  if (argc < 2) {
    std::cerr << "Missing input: number of particles to simulate." << std::endl;
    return -1;
  }
  sscanf(argv[1],"%d",&args.N);
  for (int i=2;i<argc;i++) {
    std::string token(argv[i]);
    if (token == "-CPU")
        args.device = CPU;
    else if (token == "-GPU")
        args.device = GPU;
    else if (token == "-debug")
      args.debug = true;
    else if (token == "-render")
      args.render = true;
    else {
      std::cout << "Unrecognized argument: " << token << std::endl;
       return -1;
    }
  } // Finished parsing input arguments

  std::cout << "Simulator starting." << std::endl;

  // Create Simulator object
  Simulator *sim = new Simulator(args);

  // Populate
  clock_t timer = clock();
  sim->generateParticles();
  timer = clock() - timer;
  std::cout << "Populated in " << (float)timer/CLOCKS_PER_SEC << " seconds." << std::endl;

  // Propagate
  timer = clock();
  sim->propagate();
  timer = clock() - timer;
  std::cout << "Propagated in " << (float)timer/CLOCKS_PER_SEC << " seconds." << std::endl;

  delete sim;

  std::cout << "Simulator exiting." << std::endl;
  return 0;
}