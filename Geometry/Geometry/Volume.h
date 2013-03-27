#ifndef NA63_GEOMETRY_VOLUME_H
#define NA63_GEOMETRY_VOLUME_H

#include "Material.h"

namespace na63 {

  typedef struct {
    float x, y, z;
  } ThreeVector;

  // Pointer to the struct containing parameters of the volume
  typedef void* VolumePars;

  /**
   * Abstract class. Derived classes must override Inside() and ParameterSize()
   */
  class Volume {

  public:
    Volume(Material *m, bool defined_before) {
      material_ = m;
      volume_types += defined_before;
    }
    ~Volume();

    virtual bool Inside(ThreeVector point) =0;
    virtual unsigned ParameterSize() =0;
    virtual VolumePars pars() =0;

  protected:
    Material *material() { return material_; };

  private:
    static unsigned volume_types;
    Material *material_;

  };

  unsigned Volume::volume_types = 0;

}

#endif /* NA63_GEOMETRY_VOLUME_H */