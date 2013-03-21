#ifndef NA63_SIMULATION_CUDA_HELPER_H
#define NA63_SIMULATION_CUDA_HELPER_H

#include <cuda.h>
#include <cuda_runtime.h>

namespace NA63 {

int error(cudaError_t err);

}

#endif