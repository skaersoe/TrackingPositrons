#ifndef NA63_GEOMETRY_VOLUME_H
#define NA63_GEOMETRY_VOLUME_H

#include <cstring>
#include <vector>
#include "Geometry/Library.hh"
#include "Geometry/Material.hh"
#include "Geometry/Geometry.hh"

namespace na63 {

/**
 * Abstract class. Derived classes must override Inside()
 */
class Volume {

public:
  Material *material;

  Volume(const char* mat_name, VolumeType vol_type)
      : material_name_(mat_name) {
    volume_type_ = vol_type;
    pars_.function_index = -1;
    pars_.material_index = -1;
  }
  ~Volume() {}

  std::string material_name() { return material_name_; }

  VolumePars GPU() {
    SetSpecificParameters((void*)pars_.specific);
    if (pars_.function_index < 0) throw "Missing function index";
    if (pars_.material_index < 0) throw "Missing material index";
    return pars_;
  }

  virtual bool Inside(const FourVector& position) const =0;

protected:
  virtual void SetSpecificParameters(void* parameters) =0;
  friend class Geometry;
  VolumeType volume_type() { return volume_type_; }
  void SetIndices(int mat_idx, int fun_idx) {
    pars_.material_index = mat_idx;
    pars_.function_index = fun_idx;
  }
  virtual InsideFunction inside_function() const =0;

private:
  std::string material_name_;
  VolumeType volume_type_;
  VolumePars pars_;

};

} // End namespace na63

#endif /* NA63_GEOMETRY_VOLUME_H */