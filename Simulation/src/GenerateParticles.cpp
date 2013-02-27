#include <stdio.h>

#include "HepMC/GenParticle.h"

using namespace HepMC;

__global__ void populate() {
  return;
}

void GenerateParticles(GenParticle* p, long N) {
  dim3 threadsPerBlock(256,1,1);
  dim3 blocksPerGrid((N-1)/threadsPerBlock.x+1,1,1);
  std::cout << blocksPerGrid.x << " blocks of " << threadsPerBlock.x << " threads = " << threadsPerBlock.x * blocksPerGrid.x << " total threads. " << std::endl;
}