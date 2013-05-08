#ifndef NA63_SIMULATION_TRACK_H
#define NA63_SIMULATION_TRACK_H

#include <iostream>
#include <cmath>
#include "Geometry/Library.hh"
#include "Simulation/Particle.hh"

namespace na63 {

/**
 * Track data struct, should always be aligned to 32 bytes to take advantage
 * of memory bursts.
 */
typedef struct {
  // int particle_index;
  int particle_id;    // According to Monte Carlo Particle Numbering Scheme
  int particle_index; // Set at propagation time for the GPU
  GPUFourVector position;
  GPUFourVector momentum;
  // Align to 64 bytes
  char padding[64 - 2*sizeof(int) - 2*sizeof(GPUFourVector)];
} GPUTrack;

class Track {

public:
  bool alive_;
  int particle_id;
  FourVector momentum;
  FourVector position;
  Particle *particle;
  Volume *volume;

  void Kill() {
    alive_ = false;
  }

  Track(int particle_id, FourVector pos, FourVector mom) 
      : position(pos), momentum(mom) {
    this->particle_id = particle_id;
    volume = nullptr;
    alive_ = true;
  }
  #ifdef RUNNING_CPP11
  Track(int particle_id)
      : Track(particle_id,FourVector(0,0,0,0),FourVector(0,0,0,0)) {}
  #endif

  bool alive() const { return alive_; }
  Float time() const { return position[3]; }
  Float energy() const { return momentum[3]; }
  Float charge() const { return particle->charge(); }
  Float mass() const { return particle->mass(); }
  Float momentum_magnitude() const {
    return sqrt(pow(energy(),2) - pow(mass(),2));
  } 
  Float beta() const {
    Float energy_squared = pow(energy(),2);
    return sqrt((energy_squared - pow(mass(),2)) / energy_squared);
  }
  Float gamma() const {
    return Gamma(energy(),mass());
  }

  /* Propagates the track by dl */
  void Step(const Float dl);

  GPUTrack GPU() const {
    GPUTrack retval;
    retval.particle_index = particle->index;
    retval.particle_id = particle_id;
    position.GPU(retval.position);
    momentum.GPU(retval.momentum);
    return retval;
  }

  Track& operator=(const GPUTrack& gpu_track);
  friend std::ostream& operator<<(std::ostream& os, const Track& t);

};

/*std::ostream& operator<<(std::ostream& os, const GPUTrack& t) {
  os << "Id " << t.particle_id << ", "
     << "(" << t.position[0]
     << "," << t.position[1]
     << "," << t.position[2]
     << "," << t.position[3] << "), "
     << "(" << t.momentum[0]
     << "," << t.momentum[1]
     << "," << t.momentum[2]
     << "," << t.momentum[3]
     << ")";
  return os;
}*/

} // End namespace na63

#endif /* NA63_SIMULATION_TRACK_H */