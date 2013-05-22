#ifndef NA63_SIMULATION_PARTICLE_H
#define NA63_SIMULATION_PARTICLE_H

#include <string>
#include <vector>
#include "Geometry/Library.hh"

namespace na63 {

  class Track; // Forward declaration
  typedef void (*ProcessFunction)(Track& track,const Material& material,
      const Float dl);

  typedef struct {
    int id;
    Float mass;
    // Align to 32 bytes
    char padding[32 - sizeof(int) - 2*sizeof(Float)];
  } ParticlePars;

  class Particle {

  public:
    int index; // For generation of GPU tracks

    Particle(const char* n, int id, Float mass) : name_(n) {
      id_ = id;
      mass_ = mass;
      index = -1;
    }
    ~Particle() {}

    int id() const { return id_; }
    Float mass() const { return mass_; }
    std::string name() const { return name_; }
    ParticlePars GPU() const {
      ParticlePars retval;
      retval.id = id_;
      retval.mass = mass_;
      return retval;
    }

    void RegisterProcess(ProcessFunction process) {
      processes.push_back(process);
    }

    void Query(Track& track,const Material& material,const Float dl) const {
      for (int i=0;i<processes.size();i++) {
        processes[i](track,material,dl);
      }
    }

  private:
    std::string name_;
    int id_;
    Float mass_;
    std::vector<ProcessFunction> processes;

  };  

}

#endif /* NA63_SIMULATION_PARTICLE_H */