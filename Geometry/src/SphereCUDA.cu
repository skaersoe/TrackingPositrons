#include "Geometry/Sphere.cuh"

namespace na63 {

  __host__ __device__ inline
  bool Sphere_Inside(ThreeVector point, const SpherePars sphere) {
    point.x -= sphere.center.x;
    point.y -= sphere.center.y;
    point.z -= sphere.center.z;
    return point.x*point.x + point.y*point.y + point.z*point.z
           < sphere.radius_squared;
  }

  __device__
  bool Sphere_InsideKernel(ThreeVector point, void* pars) {
    return Sphere_Inside(point,*(SpherePars*)pars);
  }

  __host__
  bool Sphere_InsideWrapper(ThreeVector point, const SpherePars sphere) {
    return Sphere_Inside(point,sphere);
  }

  InsideFunction Sphere::inside_function_ = Sphere_InsideKernel;

}