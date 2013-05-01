#ifndef NA63_SIMULATION_PARTICLE_H
#define NA63_SIMULATION_PARTICLE_H

#include <string>
#include "Geometry/Library.hh"

namespace na63 {

  typedef struct {
    int id;
    Float charge;
    Float mass;
    // Align to 32 bytes
    char padding[20];
  } ParticlePars;

  class Particle {

  public:
    Particle(const char* n, int id, Float charge, Float mass) : name_(n) {
      id_ = id;
      charge_ = charge;
      mass_ = mass;
    }
    ~Particle() {}

    int id() const { return id_; }
    Float mass() const { return mass_; }
    Float charge() const { return charge_; }
    std::string name() const { return name_; }
    ParticlePars GPU() const {
      ParticlePars retval;
      retval.id = id_;
      retval.charge = charge_;
      retval.mass = mass_;
      return retval;
    }

  private:
    std::string name_;
    int id_;
    Float mass_;
    Float charge_;

  };  

}

#endif /* NA63_SIMULATION_PARTICLE_H */