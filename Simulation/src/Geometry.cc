#include "Geometry/Geometry.hh"
#include "Geometry/Sphere.hh"

namespace na63 {

Geometry::Geometry(void) {
  volumes.push_back(nullptr);
  bounds = &volumes[0];
  material_arr_    = nullptr;
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

void Geometry::AddVolumeGeneric(Volume *volume) {
  // Material should already have been added
  int material_index = GetMaterialIndex(volume->material_name());
  if (material_index < 0) throw "Material doesn't exist.";
  int function_index = -1;
  VolumeType volume_type = volume->volume_type();
  // See if volume type exists
  // for (int i=0;i<volume_types.size();i++) {
  //   if (volume_types[i].type == volume_type) {
  //     function_index = i;
  //   }
  // }
  // if (function_index == -1) {
  //   // Add the type if it doesn't
  //   function_index = AddVolumeType(volume_type,volume->inside_function());
  // }
  function_index = volume_type;
  // Update volume parameters
  if (material_index < 0 || material_index >= materials.size()) {
    throw "Invalid index generated for volume material.";
  }
  if (function_index < 0/* || function_index >= volume_types.size()*/) {
    throw "Invalid index generated for volume type.";
  }
  volume->material = &materials[material_index];
  volume->SetIndices(material_index,function_index);
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
  GenerateVolumeArray();
  GenerateVolumeTypeArray();
}

bool Geometry::InBounds(const Track& track) {
  return volumes[0]->Inside(track.position);
}

void Geometry::Query(Track& track) {
  if (track.volume != nullptr) {
    if (track.volume->Inside(track.position)) return;
  }
  track.volume = nullptr;
  // Index 0 is bounds, so start from 1
  for (int i=1;i<volumes.size();i++) {
    if (volumes[i]->Inside(track.position)) {
      track.volume = volumes[i];
      return;
    }
  }
}

void Geometry::GenerateMaterialArray() {
  ParameterVectorToArray(&materials,&material_arr_);
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
    volume_arr_[i] = volumes[i]->GPU();
  }
}

void Geometry::DeleteParameterArrays() {
  delete material_arr_;
  delete volume_arr_;
  delete volume_type_arr_;
}

void Geometry::PrintContent() {
  std::cout << "Content of geometry: " << std::endl;
  for (int i=0;i<materials.size();i++) {
    std::cout << "Material \"" << materials[i].name() << "\"" << std::endl;
  }
  for (int i=0;i<volume_types.size();i++) {
    std::cout << "Volume type " << i << std::endl;
  }
  for (int i=0;i<volumes.size();i++) {
    std::cout << "Volume of type " << volumes[i]->volume_type()
              << " and material " << volumes[i]->material_name() << std::endl;
  }
}

int Geometry::materials_size() const {
  return materials.size();
}

} // End namespace na63