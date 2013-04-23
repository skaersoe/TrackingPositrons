#ifndef NA63_SIMULATION_PARTICLE_H
#define NA63_SIMULATION_PARTICLE_H

#include <string>

namespace na63 {

  typedef struct {
    int id;
    float charge;
    float mass;
    // Align to 32 bytes
    char padding[20];
  } ParticlePars;

  class Particle {

  public:
    Particle(const char* n, int id, float charge, float mass) : name_(n) {
      pars_.id = id;
      pars_.charge = charge;
      pars_.mass = mass;
    }
    ~Particle() {}

    int id() const { return pars_.id; }
    float mass() const { return pars_.mass; }
    float charge() const { return pars_.charge; }
    std::string name() const { return name_; }
    ParticlePars pars() const { return pars_; }

  private:
    std::string name_;
    ParticlePars pars_;

  };  

}

#endif /* NA63_SIMULATION_PARTICLE_H */