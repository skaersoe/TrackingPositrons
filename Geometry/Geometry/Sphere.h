#ifndef NA63_GEOMETRY_SPHERE_H
#define NA63_GEOMETRY_SPHERE_H

#include "Geometry/Volume.h"

namespace na63 {

  typedef struct {
    Material *material;
    ThreeVector center;
    float radius_squared;
  } SpherePars;

  class Sphere : public Volume {

  public:
    Sphere(Material *material, ThreeVector center, float radius);
    ~Sphere();

    virtual bool Inside(ThreeVector point);
    virtual unsigned ParameterSize() { return sizeof(SpherePars); }
    virtual VolumePars pars() { return (VolumePars)&pars_; }

  private:
    SpherePars pars_;

    static bool defined_before;

  };

  bool Sphere::defined_before = false;

}

#endif /* NA63_GEOMETRY_SPHERE_H */