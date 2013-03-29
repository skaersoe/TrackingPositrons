#include <iostream>

#include "Simulation/CudaHelper.hh"

namespace na63 {

  int CudaError(cudaError_t err) {
    if (err == cudaSuccess) return 0;
    std::cerr << "CUDA Error: " << cudaGetErrorString(err) << std::endl;
    return -1;
  }

}