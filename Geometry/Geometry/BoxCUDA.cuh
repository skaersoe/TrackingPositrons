#ifndef NA63_GEOMETRY_BOXCUDA_CUH
#define NA63_GEOMETRY_BOXCUDA_CUH

#include "Geometry/Box.hh"
#include "Geometry/LibraryCUDA.cuh"

namespace na63 {

__device__ 
bool Box_Inside(const GPUFourVector& point, const void* pars);

} // End namespace na63

#endif /* NA63_GEOMETRY_BOXCUDA_CUH */