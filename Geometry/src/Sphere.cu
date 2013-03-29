#include "Geometry/Sphere.hh"

namespace na63 {

  __host__ __device__ inline
  bool Sphere_InsideKernel(ThreeVector point, SpherePars pars) {
    point.x -= pars.center.x;
    point.y -= pars.center.y;
    point.z -= pars.center.z;
    return point.x*point.x + point.y*point.y + point.z*point.z
           < pars.radius_squared;
  }

}