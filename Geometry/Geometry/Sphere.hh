#ifndef NA63_GEOMETRY_SPHERE_H
#define NA63_GEOMETRY_SPHERE_H

#include "Geometry/Volume.hh"

namespace na63 {

  typedef struct {
    ThreeVector center;
    float radius_squared;
    // Pad to volume parameter size
    char padding[(int)(
      VOLUME_PARAMETER_SIZE
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
    Sphere(const char* n, ThreeVector center, float radius);
    Sphere(const Sphere& other);
    ~Sphere() {}

    virtual bool Inside(ThreeVector point) const;

  protected:
    virtual InsideFunction inside_function() const { return inside_function_; }

  private:
    static InsideFunction inside_function_;

    SpherePars *sphere_pars;

  };

}

#endif /* NA63_GEOMETRY_SPHERE_H */