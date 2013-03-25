#ifndef NA63_GEOMETRY_GEOMETRY_H
#define NA63_GEOMETRY_GEOMETRY_H

#include <vector>

#include "Geometry/Material.h"
#include "Simulation/Particle.h"

namespace na63 {

  class Geometry {

  public:
    Geometry(void);
    ~Geometry();

    void AddMaterial(Material m) {
      materials_.push_back(m);
    }
    void AddParticle(Particle p) {
      particles_.push_back(p);
    }
    /**
     * Should be called explicitly before requesting parameter arrays, but the
     * individual arrays will be generated if requested before they have been
     * generated.
     */
    void GenerateParameterArrays();
    int materials_size() { return materials_.size(); }
    int particles_size() { return particles_.size(); }
    MaterialPars *material_arr();
    ParticlePars *particle_arr();

  private:
    std::vector<Material> materials_;
    std::vector<Particle> particles_;
    MaterialPars *material_arr_;
    ParticlePars *particle_arr_;

    void GenerateMaterialArray();
    void GenerateParticleArray();
    void DeleteParameterArrays();

  };

}

#endif /* NA63_GEOMETRY_GEOMETRY_H */