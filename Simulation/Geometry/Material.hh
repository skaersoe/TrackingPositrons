#ifndef NA63_GEOMETRY_MATERIAL_H
#define NA63_GEOMETRY_MATERIAL_H

#include <string>
#include "Geometry/Library.hh"
#include "Geometry/Constants.hh"

namespace na63 {

  typedef struct {
    Float atomic_number;
    Float density;
    Float mean_excitation_potential;
    Float radiation_length;
    Float electron_density;
    Float coulomb_correction;
    char padding[32 - 5*sizeof(Float)];
  } MaterialPars;

  class Material {

  private:
    Float ComputeCoulombCorrection() const {
      static const Float k1 = 0.0083 , k2 = 0.20206 ,k3 = 0.0020 , k4 = 0.0369;
      Float az2 = (kFineStructure*atomic_number_)
         * (kFineStructure*atomic_number_);
      Float az4 = az2 * az2;
      return (k1*az4 + k2 + 1.0/(1.0+az2)) * az2 - (k3*az4 + k4)*az4;
    }

  public:
    Material(
        const char* n,
        Float atomic_number,
        Float rho,
        Float mean_excitation_potential,
        Float rl) : name_(n) {

      atomic_number_ = atomic_number;
      mean_excitation_potential_ = mean_excitation_potential;
      radiation_length_ = rl;
      density_ = rho;
      coulomb_correction_ = ComputeCoulombCorrection();

    }

    std::string name() const { return name_; }
    Float atomic_number() const { return atomic_number_; }
    Float mean_excitation_potential() const { return mean_excitation_potential_; }
    Float radiation_length() const { return radiation_length_; }
    // Since there are currently no composite materials, number of electrons is
    // just the atomic number.
    Float electron_density() const { return atomic_number_; }
    Float coulomb_correction() const { return coulomb_correction_; }
    Float density() const { return density_; }
    MaterialPars GPU() const {
      MaterialPars retval;
      retval.atomic_number = atomic_number_;
      retval.density = density_;
      retval.mean_excitation_potential = mean_excitation_potential_;
      retval.radiation_length = radiation_length_;
      // No composite materials
      retval.electron_density = atomic_number_;
      retval.coulomb_correction = coulomb_correction_;
      return retval;
    }

  private:
    std::string name_;
    Float atomic_number_;
    Float mean_excitation_potential_;
    Float radiation_length_;
    // Float electron_density_;
    Float coulomb_correction_;
    Float density_;

  };

}

#endif /* NA63_GEOMETRY_MATERIAL_H */