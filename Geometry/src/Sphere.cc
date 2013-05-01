#include "Geometry/Sphere.hh"

#include <cmath>

namespace na63 {

Sphere::Sphere(const char* n, ThreeVector center, Float radius)
      : Volume(n,SPHERE) {
  this->center = center;
  this->radius = radius;
}

Sphere::Sphere(const Sphere& other) : Volume(other) {
  this->center = other.center;
  this->radius = other.radius;
}

bool Sphere::Inside(ThreeVector point) const {
  point -= center;
  return pow(point[0],2) +
         pow(point[1],2) +
         pow(point[2],2)
         < pow(radius,2);
}

void Sphere::SetSpecificParameters(void *parameters) {
  SpherePars* sphere = (SpherePars*)parameters;
  center.GPU(sphere->center);
  sphere->radius_squared = pow(radius,2);
}

} // End namespace na63