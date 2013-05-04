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

private:
  int particle_id;
  FourVector momentum;
  FourVector position;
  Particle *particle;
  Volume *volume;
  friend class Simulator;
  friend class Geometry;

public:
  Track(int particle_id, FourVector pos, FourVector mom) 
      : position(pos), momentum(mom) {
    this->particle_id = particle_id;
  }
  #ifdef RUNNING_CPP11
  Track(int particle_id)
      : Track(particle_id,FourVector(0,0,0,0),FourVector(0,0,0,0)) {}
  #endif

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
    return energy() / mass();
  }

  /* Propagates the track by dl */
  void Step(const Float dl) {
    ThreeVector change = SphericalToCartesian(
        dl,momentum.Theta(),momentum.Phi());
    position[0] += change[0];
    position[1] += change[1];
    position[2] += change[2];
    position[3] += dl * (energy() / (momentum_magnitude() * c));
  }

  GPUTrack GPU() const {
    GPUTrack retval;
    retval.particle_id = particle_id;
    position.GPU(retval.position);
    momentum.GPU(retval.momentum);
    return retval;
  }

  Track& operator=(const GPUTrack& gpu_track) {
    particle_id = gpu_track.particle_id;
    position = gpu_track.position;
    momentum = gpu_track.momentum;
    return *this;
  }

  friend std::ostream& operator<<(std::ostream& os, const Track& t) {
    if (t.particle != NULL) {
      os << t.particle->name();
    } else {
      os << t.particle_id;
    }
    os << ", " << t.position << ", "
       << t.momentum;
    return os;
  }

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