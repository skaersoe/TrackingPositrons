#ifndef NA63_GEOMETRY_TYPES_H
#define NA63_GEOMETRY_TYPES_H

namespace na63 {

  #define VOLUME_PARAMETER_SIZE 128 - 3*sizeof(int)

  #if __cplusplus == 201103
  #define RUNNING_CPP11
  #endif
  
  // All volumes must be defined here
  typedef enum {
    SPHERE,
    BOX
  } VolumeType;

  typedef struct {
    float x, y, z;
  } ThreeVector;

  typedef struct {
    // Generic fields
    int function_index;
    int volume_index;
    int material_index;
    // Volume-specific fields. Derived classes must pad parameters to this size
    char specific[VOLUME_PARAMETER_SIZE];
  } VolumePars;

  typedef bool (*InsideFunction)(ThreeVector,void*);

  class Geometry;
  class Volume;
  class Material;

}

#endif /* NA63_GEOMETRY_TYPES_H */