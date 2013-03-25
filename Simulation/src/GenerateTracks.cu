#include <iostream>

#include "Simulation/GenerateTracks.h"
#include "Simulation/CudaHelper.h"

namespace na63 {

  __global__
  void CUDA_GenerateTracks_Electrons(Track *t_, const int N) {
    int index = threadIdx.x + blockDim.x * blockIdx.x;
    if (index >= N) return;
    float px = -0.5e-2 + 1e-2 * (index + 1) / N;
    float py =  0.5e-2 - 1e-2 * (index + 1) / N;
    float pz =  0.1 + 9.9 * (index + 1) / N;
    Track* t = &t_[index];
    t->p[0] = px;
    t->p[1] = py;
    t->p[2] = pz;
    t->r[0] = 0;
    t->r[1] = 0;
    t->r[2] = 0;
    t->particle_id = 11; // Electron
    t->particle_index = -1;
  }

  __host__
  void GenerateTracks(Track *tracks, SimulatorPars args) {
    int threadsPerBlock = 256;
    int blocksPerGrid = (args.N - 1) / threadsPerBlock + 1;
    const int dataSize = args.N*sizeof(Track);

    Track* devicePtr = NULL;
    
    if (CudaError(cudaMalloc((void**)&devicePtr,dataSize))) return;
    if (CudaError(cudaMemcpy((void*)devicePtr,tracks,dataSize,cudaMemcpyHostToDevice))) return;

    if (args.debug) {
      std::cout << "Copied " << args.N << " instances of size " << sizeof(Track) << " bytes each, resulting in a total of " << dataSize << " bytes of data on the device." << std::endl;
      std::cout << "About to initialize " << blocksPerGrid << " blocks of " << threadsPerBlock << " each, resulting in a total of " << blocksPerGrid * threadsPerBlock << " threads." << std::endl;
    }
    
    CUDA_GenerateTracks_Electrons<<<blocksPerGrid,threadsPerBlock>>>(devicePtr,args.N);
    cudaDeviceSynchronize();

    if (CudaError(cudaMemcpy(tracks,devicePtr,dataSize,cudaMemcpyDeviceToHost))) return;
    if (CudaError(cudaFree(devicePtr))) return;
  }

}