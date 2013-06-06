#include "Geometry/Sphere.hh"

#include <cmath>

namespace na63 {

Sphere::Sphere(const char* n, ThreeVector c, Float r)
      : Volume(n,SPHERE) {
  center = c;
  radius = r;
}

bool Sphere::Inside(const FourVector& position) const {
  return pow(position[0] - center[0],2) +
         pow(position[1] - center[1],2) +
         pow(position[2] - center[2],2)
         < pow(radius,2);
}

void Sphere::SetSpecificParameters(void *parameters) {
  SpherePars* sphere = (SpherePars*)parameters;
  center.GPU(sphere->center);
  sphere->radius_squared = pow(radius,2);
}

} // End namespace na63