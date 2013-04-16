#ifndef NA63_GEOMETRY_VOLUME_H
#define NA63_GEOMETRY_VOLUME_H

#include <cstring>
#include "Geometry/Material.hh"

#define VOLUME_PARAMETER_SIZE 128 - sizeof(KernelIndex) - sizeof(int)

namespace na63 {

  // This should probably be defined elsewhere
  typedef struct {
    float x, y, z;
  } ThreeVector;

  // All volumes must be defined here
  typedef enum {
    SPHERE,
    BOX
  } KernelIndex;

  // All derived classes must pad to this size
  typedef struct {
    // Generic fields
    KernelIndex kernel_index;
    int material_index;
    // Volume-specific fields
    char specific[VOLUME_PARAMETER_SIZE];
  } VolumePars;

  typedef bool (*InsideKernel)(ThreeVector,void*);

  /**
   * Abstract class. Derived classes must override Inside()
   */
  class Volume {

  public:
    Volume(Material *material, KernelIndex kernel_index) {
      material_ = material;
      pars_.kernel_index = kernel_index;
    }
    ~Volume() {}

    VolumePars pars() { return pars_; }

    virtual bool Inside(ThreeVector point) const =0;
    /**
     * Should return the static function pointer to the kernel function of the
     * given volume type.
     */
    virtual InsideKernel inside_kernel() const =0;

  protected:
    void* SpecificParameters() {
      return (void*)&pars_.specific;
    }

  private:
    Material *material_;
    VolumePars pars_;

    Material *material() const { return material_; };

  };

}

#endif /* NA63_GEOMETRY_VOLUME_H */