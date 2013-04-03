#ifndef NA63_GEOMETRY_VOLUME_H
#define NA63_GEOMETRY_VOLUME_H

#include "Geometry/Material.hh"

/* 
   Parameter size is hardcoded as we only know sizes at compile time.
   Derived classes are expected to pad to this value with something like:

     typedef struct {
       float a[2];
       char b[3];
       char padding[
         (int)(VOLUME_PARAMETER_SIZE - 2*sizeof(float) - 3*sizeof(char))
       ];
     } SomeParameters;

   REMEMBER to cast to int, as sizeof() is unsigned and will result in an
   enormous value instead of a negative number for parameter sets that exceed
   the maximum size.
   They should of course fail at compile time, and not cause the program to
   attempt to allocate monstrous amounts of memory.
*/

#define VOLUME_PARAMETER_SIZE 128

namespace na63 {

  // This should probably be defined elsewhere...
  typedef struct {
    float x, y, z;
  } ThreeVector;

  // All derived classes must pad to this size
  typedef struct {
    char size[VOLUME_PARAMETER_SIZE];
  } VolumePars;

  typedef bool (*InsideKernel)(ThreeVector,VolumePars);

  /**
   * Abstract class. Derived classes must override Inside()
   */
  class Volume {

  public:
    Volume(Material *m) {
      material_ = m;
    }
    ~Volume() {}

    virtual bool Inside(ThreeVector point) const =0;
    /**
     * Although parameters of derived classes are of the same size, C++ does not
     * allow structs to be cast directly. Instead, dirty tricks are used by
     * casting the pointer:
     *
     *   return *((VolumePars*)&pars_);
     *
     * This result is technically undefined, so correctness is not guaranteed.
     * A union could be used, but would require each derived class' parameters
     * to be declared in the union explicitly.
     *
     * For added security, a static_assert() can be added to derived classes as
     * below:
     *
     *   static_assert(sizeof(SpherePars) == sizeof(VolumePars),
     *     "Incorrect parameter size of class derived from Volume");
     */
    virtual VolumePars pars() const =0;
    /**
     * Should return the static function pointer to the kernel function of the
     * given volume type.
     */
    virtual InsideKernel kernel() const =0;

  private:
    Material *material_;

    Material *material() const { return material_; };

  };

}

#endif /* NA63_GEOMETRY_VOLUME_H */