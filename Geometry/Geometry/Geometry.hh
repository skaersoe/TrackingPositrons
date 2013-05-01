#ifndef NA63_GEOMETRY_GEOMETRY_H
#define NA63_GEOMETRY_GEOMETRY_H

#include <iostream>
#include <vector>
#ifdef RUNNING_CPP11
#include <type_traits>
#endif /* RUNNING_CPP11 */

#include "Geometry/Material.hh"
#include "Simulation/Particle.hh"
#include "Simulation/Track.hh"
#include "Geometry/Volume.hh"
#include "Geometry/Library.hh"

// Allow for easy inclusion
#include "Geometry/Sphere.hh"

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
  void AddParticle(Particle p) {
    particles.push_back(p);
  }
  /**
   * Should be called explicitly before requesting parameter arrays, but the
   * individual arrays will be generated if requested before they have been
   * generated.
   */
  void GenerateParameterArrays();
  int GetParticleIndex(int id);
  Particle* GetParticle(int id);
  int materials_size() const { return materials.size(); }
  int particles_size() const { return particles.size(); }
  int volumes_size()   const { return volumes.size();   }
  MaterialPars   *material_arr();
  ParticlePars   *particle_arr();
  VolumePars     *volume_arr();
  InsideFunction *volume_type_arr();
  /** For debugging purposes. */
  void PrintContent();
  bool InBounds(Track *t);
  void Query(Track *t);

private:
  std::vector<Material> materials;
  std::vector<Particle> particles;
  // Since volume is abstract, only pointers can be maintained here
  std::vector<Volume*>  volumes;
  std::vector<VolumeTypeFunction> volume_types;
  MaterialPars   *material_arr_;
  ParticlePars   *particle_arr_;
  VolumePars     *volume_arr_;
  InsideFunction *volume_type_arr_;

  void GenerateMaterialArray();
  void GenerateParticleArray();
  void GenerateVolumeArray();
  void GenerateVolumeTypeArray();
  void DeleteParameterArrays();
  int GetVolumeIndex(VolumeType type);
  int GetMaterialIndex(std::string material);
  int AddVolumeType(VolumeType type, InsideFunction function);
  int AddVolumeGeneric(Volume *volume);

public:
  template <class VolumeType>
  void AddVolume(VolumeType volume) {
    #ifdef RUNNING_CPP11
    if (!std::is_base_of<Volume,VolumeType>::value) {
      std::cerr << "Volume type must derive from Volume." << std::endl;
      return;
    }
    #endif
    if (AddVolumeGeneric((Volume*)&volume) != 0) return;
    volumes.push_back(new VolumeType(volume));
  }

};

} // End namespace na63

#endif /* NA63_GEOMETRY_GEOMETRY_H */