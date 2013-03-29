#include "Geometry/Geometry.hh"

namespace na63 {

  Geometry::Geometry(void) {
    material_arr_ = nullptr;
    particle_arr_ = nullptr;
    volume_arr_   = nullptr;
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

  VolumePars* Geometry::volume_arr() {
    if (volume_arr_ == nullptr)
      GenerateVolumeArray();
    return volume_arr_;
  }

  void Geometry::GenerateParameterArrays() {
    DeleteParameterArrays();
    GenerateMaterialArray();
    GenerateParticleArray();
    GenerateVolumeArray();
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

  void Geometry::GenerateVolumeArray() {
    int size = volumes_.size();
    volume_arr_ = new VolumePars[size];
    for (int i=0;i<size;i++) {
      volume_arr_[i] = volumes_[i]->pars();
    }
  }

  void Geometry::DeleteParameterArrays() {
    delete material_arr_;
    delete particle_arr_;
    delete volume_arr_;
  }

}