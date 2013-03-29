#ifndef NA63_GEOMETRY_SPHERE_H
#define NA63_GEOMETRY_SPHERE_H

#include "Geometry/Volume.hh"

namespace na63 {

  typedef struct {
    Material *material;
    ThreeVector center;
    float radius_squared;
    // Pad to volume parameter size
    char padding[(int)(
      VOLUME_PARAMETER_SIZE
      - sizeof(Material*)
      - sizeof(ThreeVector)
      - sizeof(float)
    )];
  } SpherePars;

  /**
   * Very simple geometric object to test inheritance and parameter
   * functionality. Might be nice to keep.
   */
  class Sphere : public Volume {

  public:
    Sphere(Material *material, ThreeVector center, float radius);
    ~Sphere();

    virtual bool Inside(ThreeVector point);
    virtual VolumePars pars() {
      // Dirty, dirty tricks
      return *((VolumePars*)&pars_);
    }

  private:
    SpherePars pars_;

    static bool defined_before;

  };

  bool Sphere::defined_before = false;

}

#endif /* NA63_GEOMETRY_SPHERE_H */