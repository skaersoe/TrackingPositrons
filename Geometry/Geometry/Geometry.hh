#ifndef NA63_GEOMETRY_GEOMETRY_H
#define NA63_GEOMETRY_GEOMETRY_H

#include <vector>

#include "Geometry/Material.hh"
#include "Geometry/Volume.hh"
#include "Simulation/Particle.hh"

namespace na63 {

  class Geometry {

  public:
    Geometry(void);
    ~Geometry(void);

    // For now, assume there's no reason to remove instances
    void AddMaterial(Material m) {
      materials_.push_back(m);
    }
    void AddParticle(Particle p) {
      particles_.push_back(p);
    }
    void AddVolume(Volume *v) {
      volumes_.push_back(v);
    }
    /**
     * Should be called explicitly before requesting parameter arrays, but the
     * individual arrays will be generated if requested before they have been
     * generated.
     */
    void GenerateParameterArrays();
    int materials_size() const { return materials_.size(); }
    int particles_size() const { return particles_.size(); }
    int volumes_size()   const { return volumes_.size(); }
    MaterialPars *material_arr();
    ParticlePars *particle_arr();
    VolumePars   *volume_arr();

  private:
    std::vector<Material> materials_;
    std::vector<Particle> particles_;
    // Since volume is abstract, we can only maintain pointers here
    std::vector<Volume*>  volumes_;
    MaterialPars *material_arr_;
    ParticlePars *particle_arr_;
    VolumePars   *volume_arr_;

    void GenerateMaterialArray();
    void GenerateParticleArray();
    void GenerateVolumeArray();
    void DeleteParameterArrays();

  };

}

#endif /* NA63_GEOMETRY_GEOMETRY_H */