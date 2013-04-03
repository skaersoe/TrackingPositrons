#include "Geometry/Sphere.cuh"

namespace na63 {

  __host__ __device__ inline
  bool Sphere_Inside(ThreeVector point, const SpherePars pars) {
    point.x -= pars.center.x;
    point.y -= pars.center.y;
    point.z -= pars.center.z;
    return point.x*point.x + point.y*point.y + point.z*point.z
           < pars.radius_squared;
  }

  __device__
  bool Sphere_InsideKernel(ThreeVector point, const VolumePars p) {
    return Sphere_Inside(point,*(SpherePars*)&p);
  }

  __host__
  bool Sphere_InsideWrapper(ThreeVector point, const SpherePars pars) {
    return Sphere_Inside(point,pars);
  }

  InsideKernel Sphere::inside_kernel_ = Sphere_InsideKernel;

}