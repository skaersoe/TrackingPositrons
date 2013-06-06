#ifndef NA63_GEOMETRY_SPHERECUDA_CUH
#define NA63_GEOMETRY_SPHERECUDA_CUH

#include "Geometry/Sphere.hh"
#include "Geometry/LibraryCUDA.cuh"

namespace na63 {

__device__
bool Sphere_Inside(const GPUFourVector& point, const void* parameters);

} // End namespace na63

#endif /* NA63_GEOMETRY_SPHERECUDA_CUH */