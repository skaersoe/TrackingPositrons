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
  int alive;
  int particle_id;    // According to Monte Carlo Particle Numbering Scheme
  int particle_index; // Set at propagation time for the GPU
  int volume_index;   // Remembers which volume the particle is currently inside
  int initial_index;
  GPUFourVector position;
  GPUFourVector momentum;
  // Align to 64 bytes
  char padding[64 - 5*sizeof(int) - 2*sizeof(GPUFourVector)];
} GPUTrack;

#ifdef RUNNING_CPP11
static_assert(sizeof(GPUTrack) == 64,"Unaligned GPUTrack struct");
#endif /* RUNNING_CPP11 */

class Simulator;

class Track {

public:
  bool alive;
  int particle_id;
  int initial_index;
  FourVector momentum;
  FourVector position;
  Simulator *simulator;
  Particle *particle;
  Volume *volume;
  Track *mother;

  void Kill() {
    alive = false;
  }

  Track(const int particle_id, const FourVector pos, const FourVector mom) 
      : position(pos), momentum(mom) {
    this->particle_id = particle_id;
    volume = nullptr;
    mother = nullptr;
    alive = true;
  }
  #ifdef RUNNING_CPP11
  Track(const int particle_id)
      : Track(particle_id,FourVector(0,0,0,0),FourVector(0,0,0,0)) {}
  #endif /* RUNNING_CPP11 */

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

  void SpawnChild(Track child);
  void Boost(const Float bx, const Float by, const Float bz);
  /* Propagates the track by dl */
  void Step(const Float dl);

  void UpdateMomentum(const Float change) {
    Float length = CartesianToSpherical_R(momentum[0],momentum[1],momentum[2]);
    momentum[0] += change * momentum[0] / length;
    momentum[1] += change * momentum[1] / length;
    momentum[2] += change * momentum[2] / length;
    momentum[3] += change;
  }

  GPUTrack GPU() const {
    GPUTrack retval;
    retval.alive = (alive) ? 1 : 0;
    if (particle == nullptr || particle->index < 0) {
      throw "Particle not properly registered to track.";
    }
    retval.particle_index = particle->index;
    retval.particle_id = particle_id;
    retval.volume_index = -1;
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