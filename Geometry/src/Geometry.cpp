#include "Geometry/Geometry.h"

namespace na63 {

  Geometry::Geometry(void) {
    material_arr_ = nullptr;
    particle_arr_ = nullptr;
  }

  Geometry::~Geometry() {
    DeleteParameterArrays();
  }

  MaterialPars* Geometry::material_arr() {
    if (material_arr_ == nullptr)
      GenerateMaterialArray();
    return material_arr_;
  }

  ParticlePars* Geometry::particle_arr() {
    if (particle_arr_ == nullptr)
      GenerateParticleArray();
    return particle_arr_;
  }

  void Geometry::GenerateParameterArrays() {
    DeleteParameterArrays();
    GenerateMaterialArray();
    GenerateParticleArray();
  }

  void Geometry::GenerateMaterialArray() {
    int size = materials_.size();
    material_arr_ = new MaterialPars[size];
    for (int i=0;i<size;i++) {
      material_arr_[i] = materials_[i].pars();
    }
  }

  void Geometry::GenerateParticleArray() {
    int size = particles_.size();
    particle_arr_ = new ParticlePars[size];
    for (int i=0;i<size;i++) {
      particle_arr_[i] = particles_[i].pars();
    }
  }

  void Geometry::DeleteParameterArrays() {
    delete material_arr_;
    delete particle_arr_;
  }

}