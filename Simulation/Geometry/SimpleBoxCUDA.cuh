#ifndef NA63_GEOMETRY_SIMPLEBOXCUDA_CUH
#define NA63_GEOMETRY_SIMPLEBOXCUDA_CUH

#include "Geometry/SimpleBox.hh"
#include "Geometry/LibraryCUDA.cuh"

namespace na63 {

__device__
bool SimpleBox_Inside(const GPUFourVector& point, const void* parameters);

} // End namespace na63

#endif /* NA63_GEOMETRY_SIMPLEBOXCUDA_CUH */