#ifndef NA63_GEOMETRY_SIMPLEBOX_H
#define NA63_GEOMETRY_SIMPLEBOX_H

#include "Geometry/Volume.hh"

namespace na63 {

typedef struct {
  GPUThreeVector center;
  GPUThreeVector size;
  // Pad to volume parameter size
  char padding[(int)(
    VOLUME_PARAMETER_SIZE
    - 2*sizeof(GPUThreeVector)
  )];
} SimpleBoxPars;

#ifdef RUNNING_CPP11
// Some extra safety to shield against human errors
static_assert(sizeof(SimpleBoxPars) == VOLUME_PARAMETER_SIZE,
    "Incorrect parameter size of class derived from Volume");
#endif

class SimpleBox : public Volume {

public:
  SimpleBox(const char* n, ThreeVector c, ThreeVector s)
      : Volume(n,SIMPLEBOX), center(c), size(s) {}

  virtual bool Inside(const FourVector& position) const;

protected:
  virtual void SetSpecificParameters(void *parameters);
  virtual InsideFunction inside_function() const { return inside_function_; }

private:
  ThreeVector center;
  ThreeVector size;
  static InsideFunction inside_function_;

};

} // End namespace na63

#endif /* NA63_GEOMETRY_SIMPLEBOX_H */