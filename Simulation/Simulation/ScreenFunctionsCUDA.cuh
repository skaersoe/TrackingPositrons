#ifndef NA63_SIMULATION_SCREENFUNCTIONSCUDA_CUH
#define NA63_SIMULATION_SCREENFUNCTIONSCUDA_CUH

#include "Geometry/LibraryCUDA.cuh"

namespace na63 {

__device__ Float CUDA_ScreenFunction1(Float screen_variable);
__device__ Float CUDA_ScreenFunction2(Float screen_variable);

} // End namespace na63

#endif /* NA63_SIMULATION_SCREENFUNCTIONSCUDA_CUH */