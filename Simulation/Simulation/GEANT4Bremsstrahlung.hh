#ifndef NA63_SIMULATION_GEANT4BREMSSTRAHLUNG_H
#define NA63_SIMULATION_GEANT4BREMSSTRAHLUNG_H

#include <cfloat>

#include "Geometry/Library.hh"
#include "Simulation/Track.hh"
#include "Simulation/Simulator.hh"
#include "Simulation/Particle.hh"
#include "Simulation/Process.hh"

namespace na63 {

Float ComputeParametrizedDXSectionPerAtom(Float kinetic_energy, Float gamma_energy, Float Z);

class Bremsstrahlung : public Process {

public:

  Bremsstrahlung(const Particle* p, const Float threshold = 2.0 * kElectronMass);

  void Query(Track* track, const Material* material,
      const Float dl);

  void SetMaterial(const Material* material_new);
  void SetParticle(const Particle* particle_new);
  void SetSecondaryThreshold(const Float threshold);
  void SetupForMaterial(const Material* mat, const Float kin_energy);

  Float ComputeDXSectionPerAtom(const Float gamma_energy) const;
  Float ComputeXSectionPerAtom(const Float cut) const;
  Float ComputeCrossSectionPerAtom(const Float cut_energy,
      const Float max_energy = DBL_MAX) const;
  void SampleSecondaries(
    const Material* couple,
    Track* track,
    Float cut_energy = 2*kElectronMass,
    Float max_energy = FLT_MAX
  );

private:

  Float f_max;
  Float Z, lnZ, z13, z23;
  Float fel, finel;
  Float coulomb_correction;
  const Particle *particle;
  Float particle_mass;
  const Material *material;
  bool is_electron;
  Float density_factor;
  Float kinetic_energy;
  Float total_energy;
  Float density_correction;
  Float min_threshold;
  Float low_kinetic_energy;
  Float low_limit;
  Float secondary_threshold;

};

const Float xgi[] = { 0.0199, 0.1017, 0.2372, 0.4083, 0.5917,
    0.7628, 0.8983, 0.9801 };
const Float wgi[] = { 0.0506, 0.1112, 0.1569, 0.1813, 0.1813,
    0.1569, 0.1112, 0.0506 };

} // End namespace na63

#endif /* NA63_SIMULATION_GEANT4BREMSSTRAHLUNG_H */