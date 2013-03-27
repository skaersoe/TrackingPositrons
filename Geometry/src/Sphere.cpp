#include <cmath>

#include "Geometry/Sphere.h"

namespace na63 {

  Sphere::Sphere(Material *material, ThreeVector center, float radius)
        : Volume(material, defined_before) {
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