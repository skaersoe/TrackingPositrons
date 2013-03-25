#ifndef NA63_SIMULATION_CUDA_HELPER_H
#define NA63_SIMULATION_CUDA_HELPER_H

#include <cuda.h>
#include <cuda_runtime.h>

namespace na63 {

  int CudaError(cudaError_t err);

}

#endif