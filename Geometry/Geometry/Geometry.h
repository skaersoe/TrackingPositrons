#ifndef NA63_GEOMETRY_H
#define NA63_GEOMETRY_H

#include <vector>

#include "Geometry/Material.h"
#include "Simulation/Particle.h"

namespace na63 {

  class Geometry {

    Geometry(void);

    public:
      void AddMaterial(Material m) {
        materials.push_back(m);
      }
      void AddParticle(Particle p) {
        particles.push_back(p);
      }

    private:
      std::vector<Material> materials;
      std::vector<Particle> particles;

  };

}

#endif /* NA63_GEOMETRY_H */