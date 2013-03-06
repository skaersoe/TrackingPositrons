#include <iostream>

#include "Simulation/CudaHelper.h"

int error(cudaError_t err) {
  if (err == cudaSuccess) return 0;
  std::cerr << "CUDA Error: " << cudaGetErrorString(err) << std::endl;
  return -1;
}