#include <stdio.h>

#include "Simulation/GenerateParticles.h"

__global__ void cudaPopulateElectrons(simple_particle_t* _p, const int N) {
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index > N) return;
  float px = -0.05 + 0.1 * (N - index) / index;
  float py =  0.05 - 0.1 * (N - index) / index;
  float pz =  0.10 + 9.9 * (N - index) / index;
  simple_particle_t *p = &_p[index];
  p->p[0] = px;
  p->p[1] = py;
  p->p[2] = pz;
  p->r[0] = 0;
  p->r[1] = 0;
  p->r[2] = 0;
  p->m = 1;
  p->q = -1;
}

int error(cudaError_t err) {
  if (err == cudaSuccess) return 0;
  printf("An error occurred.\n");
  return -1;
}

void GenerateParticles(simple_particle_t* p, const int N) {
  int threadsPerBlock = 256;
  int blocksPerGrid = (N - 1) / threadsPerBlock + 1;
  const int dataSize = N*sizeof(simple_particle_t);

  void* devicePtr = (void*)p;
  
  if (error(cudaMalloc(&devicePtr,dataSize))) return;
  if (error(cudaMemcpy(devicePtr,p,dataSize,cudaMemcpyHostToDevice))) return;

  cudaPopulateElectrons<<<blocksPerGrid,threadsPerBlock>>>(p,N);

  if (error(cudaMemcpy(p,devicePtr,dataSize,cudaMemcpyDeviceToHost))) return;
  if (error(cudaFree(devicePtr))) return;
}