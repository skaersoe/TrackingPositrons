#include <iostream>

// External
#include "HepMC/GenParticle.h"

// Local
#include "GenerateParticles.h"

#define SIMULATOR_PARTICLES 1024

int main(int argc,char *argv[]) {
  std::cout << "Simulator starting." << std::endl;

  const long N = SIMULATOR_PARTICLES;
  // Allocate particles
  HepMC::GenParticle p[N];
  // Populate on GPU
  GenerateParticles(p,N);

  std::cout << "Simulator exiting." << std::endl;
  return 0;
}