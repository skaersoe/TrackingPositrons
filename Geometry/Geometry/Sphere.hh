#ifndef NA63_GEOMETRY_SPHERE_H
#define NA63_GEOMETRY_SPHERE_H

#include "Geometry/Volume.hh"

namespace na63 {

typedef struct {
  GPUThreeVector center;
  Float radius_squared;
  // Pad to volume parameter size
  char padding[(int)(
    VOLUME_PARAMETER_SIZE
    - sizeof(GPUThreeVector)
    - sizeof(Float)
  )];
} SpherePars;

#ifdef RUNNING_CPP11
// Some extra safety to shield against human errors
static_assert(sizeof(SpherePars) == VOLUME_PARAMETER_SIZE,
    "Incorrect parameter size of class derived from Volume");
#endif

/**
 * Very simple geometric object to test inheritance and parameter
 * functionality. Might be nice to keep.
 */
class Sphere : public Volume {

public:
  Sphere(const char* n, ThreeVector center, Float radius);
  Sphere(const Sphere& other);
  ~Sphere() {}

  virtual bool Inside(const FourVector& position) const;

protected:
  virtual void SetSpecificParameters(void *parameters);
  virtual InsideFunction inside_function() const { return inside_function_; }

private:
  ThreeVector center;
  Float radius;
  static InsideFunction inside_function_;

};

} // End namespace na63

#endif /* NA63_GEOMETRY_SPHERE_H */