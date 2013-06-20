#ifndef NA63_SIMULATION_PARTICLE_H
#define NA63_SIMULATION_PARTICLE_H

#include <string>
#include <vector>
#include "Geometry/Library.hh"
#include "Simulation/Process.hh"

namespace na63 {

class Track; // Forward declaration

typedef struct {
  int id;
  Float mass;
  // Align to 32 bytes
  char padding[32 - 1*sizeof(int) - 1*sizeof(Float)];
} ParticlePars;

#ifdef RUNNING_CPP11
static_assert(sizeof(ParticlePars) == 32,"Unaligned ParticlePars struct");
#endif /* RUNNING_CPP11 */

class Particle {

public:
  int index; // For generation of GPU tracks

  Particle(const char* n, int id, Float mass);

  int id() const { return id_; }
  Float mass() const { return mass_; }
  std::string name() const { return name_; }
  ParticlePars GPU() const;

  void RegisterProcess(Process *process);

  void Query(Track* track, const Material* material, const Float dl) const;

private:
  std::string name_;
  int id_;
  Float mass_;
  std::vector<Process*> processes;

};  

} // End namespace na63

#endif /* NA63_SIMULATION_PARTICLE_H */