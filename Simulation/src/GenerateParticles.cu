#include <stdio.h>

#include "GenerateParticles.h"

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

void GenerateParticles(simple_particle_t* p, const int N) {
  dim3 threadsPerBlock(256,1,1);
  dim3 blocksPerGrid((N-1)/threadsPerBlock.x+1,1,1);
  const int dataSize = N*sizeof(simple_particle_t);

  void* devicePtr = (void*)p;
  
  cudaMalloc(&devicePtr,dataSize);
  cudaMemcpy(devicePtr,p,dataSize,cudaMemcpyHostToDevice);

  cudaPopulateElectrons<<<blocksPerGrid,threadsPerBlock>>>(p,N);

  cudaMemcpy(p,devicePtr,dataSize,cudaMemcpyDeviceToHost);
  cudaFree(devicePtr);
}