#include <iostream>
#include <time.h>

#include "Simulation/GenerateParticles.h"
#include "Simulation/Particle.h"
#include "Simulation/Propagate.h"
#include "Simulation/LaunchArguments.h"

int main(int argc,char *argv[]) {

  /* Set some default values. N must be specified by user */
  launch_args_t args = { .device = GPU, .debug = false, .N=-1 };
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
    else {
      std::cout << "Unrecognized argument: " << token << std::endl;
       return -1;
    }
  } // Finished parsing input arguments

  std::cout << "Simulator starting." << std::endl;

  // Allocate particles
  simple_particle_t p[args.N];
  // Populate
  clock_t timer = clock();
  GenerateParticles(p,args);
  timer = clock() - timer;
  std::cout << "Populated in " << (float)timer/CLOCKS_PER_SEC << " seconds." << std::endl;

  // Propagate
  timer = clock();
  propagate(p,args);
  timer = clock() - timer;
  std::cout << "Propagated in " << (float)timer/CLOCKS_PER_SEC << " seconds." << std::endl;

  std::cout << "Simulator exiting." << std::endl;
  return 0;
}