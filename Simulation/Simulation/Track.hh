#ifndef NA63_SIMULATION_TRACK_H
#define NA63_SIMULATION_TRACK_H

#include <iostream>
#include <cmath>
#include "Geometry/Library.hh"
#include "Simulation/Particle.hh"

namespace na63 {

typedef enum {
  ALIVE = 0,
  DEAD = 1,
  FREE = 2,
  WAITING = 3
} TrackState;

/**
 * Track data struct, should always be aligned to 32 bytes to take advantage
 * of memory bursts.
 */
typedef struct {
  // int particle_index;
  TrackState state;
  int particle_id;    // According to Monte Carlo Particle Numbering Scheme
  int particle_index; // Set at propagation time for the GPU
  int volume_index;   // Remembers which volume the particle is currently inside
  int initial_index;  // For reconciliation with original indices
  int charge;
  GPUFourVector vertex;
  GPUFourVector position;
  GPUFourVector momentum;
  // Align to 128 bytes
  char padding[128 - 6*sizeof(int) - 3*sizeof(GPUFourVector)];
} GPUTrack;

#ifdef RUNNING_CPP11
static_assert(sizeof(GPUTrack) == 128,"Unaligned GPUTrack struct");
#endif /* RUNNING_CPP11 */

class Simulator;

class Track {

private:
  FourVector vertex_;
  int charge_;

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

  void Kill();

  void Stop();

  Track(const int particle_id, const int c, const FourVector pos,
      const FourVector mom) : position(pos), vertex_(pos), momentum(mom) {
    this->particle_id = particle_id;
    charge_ = c;
    volume = nullptr;
    mother = nullptr;
    alive = true;
  }
  #ifdef RUNNING_CPP11
  Track(const int particle_id, const int c)
      : Track(particle_id,c,FourVector(0,0,0,0),FourVector(0,0,0,0)) {}
  Track() : Track(0,0) {}
  #endif /* RUNNING_CPP11 */

  Float time() const { return position[3]; }
  Float energy() const { return momentum[3]; }
  Float kinetic_energy() const { return momentum[3] - mass(); }
  Float charge() const { return charge_; }
  Float mass() const { return particle->mass(); }
  Float momentum_magnitude() const {
    return sqrt(pow(energy(),2) - pow(mass(),2));
  }
  Float beta() const;
  ThreeVector beta_vector() const;
  Float gamma() const;
  Float theta() const {
    return momentum.Theta();
  }
  Float phi() const {
    return momentum.Phi();
  }
  FourVector vertex() const {
    return vertex_;
  }

  void SpawnChild(Track child);
  void SpawnChildren(std::vector<Track> children);
  void Boost(const Float bx, const Float by, const Float bz);
  #ifdef RUNNING_CPP11
  void Boost(ThreeVector tv) {
    Boost(tv[0],tv[1],tv[2]);
  }
  #endif
  /* Propagates the track by dl */
  void Step(const Float dl);

  void UpdateEnergy(const Float change);
  void UpdateMomentum(const FourVector change);
  void SetMomentum(const FourVector momentum_new);

  GPUTrack GPU() const {
    GPUTrack retval;
    retval.state = (alive) ? ALIVE : DEAD;
    if (particle == nullptr || particle->index < 0) {
      throw "Particle not properly registered to track.";
    }
    retval.particle_index = particle->index;
    retval.charge = charge_;
    retval.particle_id = particle_id;
    retval.volume_index = -1;
    position.GPU(retval.position);
    momentum.GPU(retval.momentum);
    vertex_.GPU(retval.vertex);
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