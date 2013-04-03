#ifndef NA63_SIMULATION_PARTICLE_H
#define NA63_SIMULATION_PARTICLE_H

namespace na63 {

  typedef struct {
    int id; // According to the Monte Carlo Particle Numbering Scheme
    int charge;
    float mass;
    // Align to 32 bytes
    char padding[24];
  } ParticlePars;

  class Particle {

  public:
    Particle(int id, float charge, float mass) {
      pars_.id = id;
      pars_.charge = charge;
      pars_.mass = mass;
    }
    ~Particle() {}

    int id() const { return pars_.id; }
    float mass() const { return pars_.mass; }
    float charge() const { return pars_.charge; }
    ParticlePars pars() const { return pars_; }

  private:
    ParticlePars pars_;

  };  

}

#endif /* NA63_SIMULATION_PARTICLE_H */