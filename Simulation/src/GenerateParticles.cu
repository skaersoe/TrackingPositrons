#include <iostream>

#include "Simulation/GenerateParticles.h"
#include "Simulation/CudaHelper.h"

namespace Simulation {

  __global__
  void cuda_generateParticles_electrons(simple_particle_t* p_, const int N) {
    int index = threadIdx.x + blockDim.x * blockIdx.x;
    if (index >= N) return;
    float px = -0.5e-2 + 1e-2 * (index + 1) / N;
    float py =  0.5e-2 - 1e-2 * (index + 1) / N;
    float pz =  0.1 + 9.9 * (index + 1) / N;
    simple_particle_t* p = &p_[index];
    p->p[0] = px;
    p->p[1] = py;
    p->p[2] = pz;
    p->r[0] = 0;
    p->r[1] = 0;
    p->r[2] = 0;
    p->m = 5.109989e-4;
    p->q = -1;
  }

  __host__
  void generateParticles(simple_particle_t* p, simulator_args_t args) {
    const int N = args.nParticles;
    int threadsPerBlock = 256;
    int blocksPerGrid = (N - 1) / threadsPerBlock + 1;
    const int dataSize = N*sizeof(simple_particle_t);

    simple_particle_t* devicePtr = NULL;
    
    if (error(cudaMalloc((void**)&devicePtr,dataSize))) return;
    if (error(cudaMemcpy((void*)devicePtr,p,dataSize,cudaMemcpyHostToDevice))) return;

    if (args.debug) {
      std::cout << "Copied " << N << " instances of size " << sizeof(simple_particle_t) << " bytes each, resulting in a total of " << dataSize << " bytes of data on the device." << std::endl;
      std::cout << "About to initialize " << blocksPerGrid << " blocks of " << threadsPerBlock << " each, resulting in a total of " << blocksPerGrid * threadsPerBlock << " threads." << std::endl;
    }
    
    cuda_generateParticles_electrons<<<blocksPerGrid,threadsPerBlock>>>(devicePtr,N);
    cudaDeviceSynchronize();

    if (error(cudaMemcpy(p,devicePtr,dataSize,cudaMemcpyDeviceToHost))) return;
    if (error(cudaFree(devicePtr))) return;
  }

}