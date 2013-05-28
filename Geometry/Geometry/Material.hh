#ifndef NA63_GEOMETRY_MATERIAL_H
#define NA63_GEOMETRY_MATERIAL_H

#include <string>
#include "Geometry/Library.hh"

namespace na63 {

  typedef struct {
    Float atomic_number;
    Float mean_excitation_potential;
    Float radiation_length;
    Float electron_density;
    Float coulomb_correction;
    char padding[32 - 5*sizeof(Float)];
  } MaterialPars;

  class Material {

  public:
    Material(const char* n, Float atomic_number,
        Float mean_excitation_potential, Float rl, Float elec_dens,
        Float coulomb_corr) : name_(n) {
      atomic_number_ = atomic_number;
      mean_excitation_potential_ = mean_excitation_potential;
      radiation_length_ = rl;
      electron_density_ = elec_dens;
      coulomb_correction_ = coulomb_corr;
    }
    ~Material() {}

    std::string name() const { return name_; }
    Float atomic_number() const { return atomic_number_; }
    Float mean_excitation_potential() const { return mean_excitation_potential_; }
    Float radiation_length() const { return radiation_length_; }
    Float electron_density() const { return electron_density_; }
    Float coulomb_correction() const { return coulomb_correction_; }
    MaterialPars GPU() const {
      MaterialPars retval;
      retval.atomic_number = atomic_number_;
      retval.mean_excitation_potential = mean_excitation_potential_;
      retval.radiation_length = radiation_length_;
      retval.electron_density = electron_density_;
      retval.coulomb_correction = coulomb_correction_;
      return retval;
    }

  private:
    std::string name_;
    Float atomic_number_;
    Float mean_excitation_potential_;
    Float radiation_length_;
    Float electron_density_;
    Float coulomb_correction_;

  };

}

#endif /* NA63_GEOMETRY_MATERIAL_H */