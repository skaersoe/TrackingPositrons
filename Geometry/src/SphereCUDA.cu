#include "Geometry/Sphere.hh"

namespace na63 {

__device__
bool Sphere_Inside(GPUThreeVector point, void* parameters) {
  SpherePars *sphere = (SpherePars*)parameters;
  point[0] -= sphere->center[0];
  point[1] -= sphere->center[1];
  point[2] -= sphere->center[2];
  return point[0]*point[0] + point[1]*point[1] + point[2]*point[2]
         < sphere->radius_squared;
}

InsideFunction Sphere::inside_function_ = Sphere_Inside;

} // End namespace na63