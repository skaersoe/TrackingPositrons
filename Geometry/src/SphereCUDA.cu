#include "Geometry/LibraryCUDA.cuh"
#include "Geometry/Sphere.hh"

namespace na63 {

__device__
bool Sphere_Inside(const GPUFourVector& point, const void* parameters) {
  SpherePars *sphere = (SpherePars*)parameters;
  return pow(point[0] - sphere->center[0],2) +
         pow(point[1] - sphere->center[1],2) +
         pow(point[2] - sphere->center[2],2)
         < sphere->radius_squared;
}

InsideFunction Sphere::inside_function_ = Sphere_Inside;

} // End namespace na63