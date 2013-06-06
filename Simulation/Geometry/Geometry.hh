#ifndef NA63_GEOMETRY_GEOMETRY_H
#define NA63_GEOMETRY_GEOMETRY_H

#include <iostream>
#include <vector>
#ifdef RUNNING_CPP11
#include <type_traits>
#endif /* RUNNING_CPP11 */

#include "Geometry/Material.hh"
#include "Simulation/Track.hh"
#include "Geometry/Volume.hh"
#include "Geometry/Library.hh"

/*
// Allow for easy inclusion
#include "Geometry/Sphere.hh"
#include "Geometry/Box.hh"
*/

namespace na63 {

typedef struct {
  VolumeType type;
  InsideFunction function;
} VolumeTypeFunction;

class Geometry {

public:
  Geometry(void);
  ~Geometry(void);

  // For now, assume there's no reason to remove instances
  void AddMaterial(Material m) {
    materials.push_back(m);
  }
  /**
   * Should be called explicitly before requesting parameter arrays, but the
   * individual arrays will be generated if requested before they have been
   * generated.
   */
  void GenerateParameterArrays();
  int materials_size() const;
  int volumes_size()   const { return volumes.size();   }
  int volume_types_size() const { return volume_types.size(); }
  MaterialPars   *material_arr();
  VolumePars     *volume_arr();
  InsideFunction *volume_type_arr();
  /** For debugging purposes. */
  void PrintContent();
  bool InBounds(const Track& t);
  void Query(Track& t);

private:
  Volume** bounds;
  std::vector<Material> materials;
  // Since volume is abstract, only pointers can be maintained here
  std::vector<Volume*>  volumes;
  std::vector<VolumeTypeFunction> volume_types;
  MaterialPars   *material_arr_;
  VolumePars     *volume_arr_;
  InsideFunction *volume_type_arr_;

  void GenerateMaterialArray();
  void GenerateVolumeArray();
  void GenerateVolumeTypeArray();
  void DeleteParameterArrays();
  int GetVolumeIndex(VolumeType type);
  int GetMaterialIndex(std::string material);
  int AddVolumeType(VolumeType type, InsideFunction function);
  void AddVolumeGeneric(Volume *volume);

public:
  template <class VolumeType>
  void AddVolume(VolumeType volume) {
    #ifdef RUNNING_CPP11
    if (!std::is_base_of<Volume,VolumeType>::value) {
      std::cerr << "Volume type must derive from Volume." << std::endl;
      return;
    }
    #endif
    AddVolumeGeneric((Volume*)&volume);
    volumes.push_back(new VolumeType(volume));
  }
  template <class VolumeType>
  void SetBounds(VolumeType volume) {
    #ifdef RUNNING_CPP11
    if (!std::is_base_of<Volume,VolumeType>::value) {
      std::cerr << "Volume type must derive from Volume." << std::endl;
      return;
    }
    #endif
    AddVolumeGeneric((Volume*)&volume);
    *bounds = new VolumeType(volume);
  }

};

} // End namespace na63

#endif /* NA63_GEOMETRY_GEOMETRY_H */