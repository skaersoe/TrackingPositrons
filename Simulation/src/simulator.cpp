#include <iostream>
#include <time.h>

// #include "HepMC/GenParticle.h"

#include "Simulation/GenerateParticles.h"
#include "Simulation/Particle.h"
#include "Simulation/Step.h"

#define SIMULATOR_PARTICLES 1024

int main(int argc,char *argv[]) {
  std::cout << "Simulator starting." << std::endl;


  const int N = SIMULATOR_PARTICLES;
  // Allocate particles
  simple_particle_t p[N];
  // Populate
  clock_t timer = clock();
  GenerateParticles(p,N);
  timer = clock() - timer;
  std::cout << "Populated in " << (float)timer/CLOCKS_PER_SEC << " seconds." << std::endl;

  std::cout << "Simulator exiting." << std::endl;
  return 0;
}