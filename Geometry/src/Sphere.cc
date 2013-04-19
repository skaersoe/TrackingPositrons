#include "Geometry/Sphere.hh"
#include "Geometry/Sphere.cuh"

namespace na63 {

  // Some extra safety to shield against human errors
  // Put here and not in header to allow header inclusion in CUDA files.
  static_assert(sizeof(SpherePars) == VOLUME_PARAMETER_SIZE,
      "Incorrect parameter size of class derived from Volume");

  Sphere::Sphere(const char* n, ThreeVector center, float radius)
        : Volume(n,SPHERE) {
    sphere_pars = (SpherePars*)SpecificParameters();
    sphere_pars->center = center;
    sphere_pars->radius_squared = radius*radius;
  }

  Sphere::Sphere(const Sphere& other) : Volume(other) {
    sphere_pars = (SpherePars*)SpecificParameters();
    *sphere_pars = *other.sphere_pars;
  }

  bool Sphere::Inside(ThreeVector point) const {
    return Sphere_InsideWrapper(point,*sphere_pars);
  }

}