#include "Geometry/Geometry.hh"

namespace na63 {

  Geometry::Geometry(void) {
    material_arr_    = nullptr;
    particle_arr_    = nullptr;
    volume_arr_      = nullptr;
    volume_type_arr_ = nullptr;
  }

  Geometry::~Geometry() {
    for (int i=0;i<volumes.size();i++) {
      delete volumes[i];
    }
    DeleteParameterArrays();
  }

  int Geometry::AddVolumeType(VolumeType type, InsideFunction function) {
    // Save the type for reconciliation. Indices will be used on the GPU
    VolumeTypeFunction type_function = {
      .type = type,
      .function = function
    };
    volume_types.push_back(type_function);
    return volume_types.size() - 1;
  }

  int Geometry::AddVolumeGeneric(Volume *volume) {
    // Material should already have been added
    int material_index = GetMaterialIndex(volume->material_name());
    if (material_index == -1) {
      std::cerr << "Material not found. Volume was not added." << std::endl;
      return -1;
    }
    int volume_index = -1;
    VolumeType volume_type = volume->volume_type();
    // See if volume type exists
    for (int i=0;i<volume_types.size();i++) {
      if (volume_types[i].type == volume_type) {
        volume_index = i;
      }
    }
    if (volume_index == -1) {
      // Add the type if it doesn't
      volume_index = AddVolumeType(volume_type,volume->inside_function());
    }
    // Update volume parameters
    volume->SetIndices(material_index,volume_index);
    return 0;
  }

  int Geometry::GetMaterialIndex(std::string material_name) {
    for (int i=0;i<materials.size();i++) {
      if (materials[i].name().compare(material_name) == 0) {
        return i;
      }
    }
    return -1;
  }

  MaterialPars* Geometry::material_arr() {
    if (material_arr_ == nullptr) {
      GenerateMaterialArray();
    }
    return material_arr_;
  }

  ParticlePars* Geometry::particle_arr() {
    if (particle_arr_ == nullptr) {
      GenerateParticleArray();
    }
    return particle_arr_;
  }

  VolumePars* Geometry::volume_arr() {
    if (volume_arr_ == nullptr) {
      GenerateVolumeArray();
    }
    return volume_arr_;
  }

  InsideFunction* Geometry::volume_type_arr() {
    if (volume_type_arr_ == nullptr) {
      GenerateVolumeTypeArray();
    }
    return volume_type_arr_;
  }

  void Geometry::GenerateParameterArrays() {
    DeleteParameterArrays();
    GenerateMaterialArray();
    GenerateParticleArray();
    GenerateVolumeArray();
    GenerateVolumeTypeArray();
  }

  template <class VectorType, class ArrayType>
  void ParameterVectorToArray(VectorType *vec, ArrayType *arr) {
    int size = vec->size();
    arr = new ArrayType[size];
    for (int i=0;i<size;i++) {
      arr[i] = (*vec)[i].pars();
    }
  }

  void Geometry::GenerateMaterialArray() {
    ParameterVectorToArray(&materials,material_arr_);
  }

  void Geometry::GenerateParticleArray() {
    ParameterVectorToArray(&particles,particle_arr_);
  }

  void Geometry::GenerateVolumeTypeArray() {
    int size = volume_types.size();
    volume_type_arr_ = new InsideFunction[size];
    for (int i=0;i<size;i++) {
      volume_type_arr_[i] = volume_types[i].function;
    }
  }

  void Geometry::GenerateVolumeArray() {
    int size = volumes.size();
    volume_arr_ = new VolumePars[size];
    for (int i=0;i<size;i++) {
      volume_arr_[i] = volumes[i]->pars();
    }
  }

  void Geometry::DeleteParameterArrays() {
    delete material_arr_;
    delete particle_arr_;
    delete volume_arr_;
    delete volume_type_arr_;
  }

  void Geometry::PrintContent() {
    std::cout << "Content of geometry: " << std::endl;
    for (int i=0;i<materials.size();i++) {
      std::cout << "Material \"" << materials[i].name() << "\"" << std::endl;
    }
    for (int i=0;i<particles.size();i++) {
      std::cout << "Particle \"" << particles[i].name() << "\"" << std::endl;
    }
    for (int i=0;i<volume_types.size();i++) {
      std::cout << "Volume type " << i << std::endl;
    }
    for (int i=0;i<volumes.size();i++) {
      std::cout << "Volume of type " << volumes[i]->volume_type()
                << " and material " << volumes[i]->material_name() << std::endl;
    }
  }

}