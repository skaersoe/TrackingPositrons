#include <cmath>

#include "Geometry/Sphere.hh"

namespace na63 {

  // Some extra safety to shield against human errors
  // Put here and not in header to allow header inclusion in CUDA files.
  static_assert(sizeof(SpherePars) == sizeof(VolumePars),
    "Incorrect parameter size of class derived from Volume");

  Sphere::Sphere(Material *material, ThreeVector center, float radius)
        : Volume(material) {
      pars_.center = center;
      pars_.radius_squared = radius*radius;
      defined_before = true;
    };

  bool Sphere::Inside(ThreeVector point) {
    point.x -= pars_.center.x;
    point.y -= pars_.center.y;
    point.z -= pars_.center.z;
    return point.x*point.x + point.y*point.y + point.z*point.z
           < pars_.radius_squared;
  }

}