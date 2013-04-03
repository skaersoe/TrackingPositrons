#ifndef NA63_SIMULATION_CUDAHELPER_H
#define NA63_SIMULATION_CUDAHELPER_H

#include <cuda.h>
#include <cuda_runtime.h>

namespace na63 {

  int CudaError(cudaError_t err);

}

#endif /* NA63_SIMULATION_CUDAHELPER_H */