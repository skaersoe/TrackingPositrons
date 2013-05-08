#ifndef NA63_GEOMETRY_MATERIAL_H
#define NA63_GEOMETRY_MATERIAL_H

#include <string>
#include "Geometry/Library.hh"

namespace na63 {

// Z_FE = 26
// m_e = 0.51098892 MeV
// I_FE = 286
  typedef struct {
    Float atomic_number;
    Float mean_excitation_potential;
    // char padding[0 - 2*sizeof(Float)];
  } MaterialPars;

  class Material {

  public:
    Material(const char* n, Float atomic_number,
        Float mean_excitation_potential) : name_(n) {
      atomic_number_ = atomic_number;
      mean_excitation_potential_ = mean_excitation_potential;
    }
    ~Material() {}

    std::string name() const { return name_; }
    Float atomic_number() const { return atomic_number_; }
    Float mean_excitation_potential() const { return mean_excitation_potential_; }
    MaterialPars GPU() const {
      MaterialPars retval;
      retval.atomic_number = atomic_number_;
      retval.mean_excitation_potential = mean_excitation_potential_;
      return retval;
    }

  private:
    std::string name_;
    Float atomic_number_;
    Float mean_excitation_potential_;

  };

}

#endif /* NA63_GEOMETRY_MATERIAL_H */