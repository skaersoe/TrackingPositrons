#ifndef NA63_SIMULATION_SIMULATOR_H
#define NA63_SIMULATION_SIMULATOR_H

#include <vector>

#include "Simulation/Track.hh"
#include "Simulation/Particle.hh"
#include "Geometry/Geometry.hh"

namespace na63 {

enum Device {CPU,GPU};

class Simulator {

public:
  // Allow these parameters to be altered freely
  Device device;
  bool debug;
  int cpu_threads;
  Float step_size;
  Geometry *geometry;

  Simulator(void);
  Simulator(Geometry *geometry);
  ~Simulator();

  int GetParticleIndex(int id) const;
  Particle* GetParticle(int id);
  void AddParticle(Particle p) {
    p.index = particles.size();
    particles.push_back(p);
  }
  int particles_size() const { return particles.size(); }
  ParticlePars *particle_arr();

  Track GetTrack(int index) const;
  void AddTracks(std::vector<Track> t);
  void ClearTracks();
  //void SetTracks(std::vector<Track> t);
  unsigned TrackSize() { return tracks.size(); }
  GPUTrack* GPUTracks();
  void CopyBackTracks();

  /**
   * Generates some hardcoded electrons.
   */
  // void GenerateTracks();
  /**
   * Runs the propagation.
   */
  void Propagate();
  void PropagateTrack(Track& track);

private:
  std::vector<Track> tracks;
  std::vector<Particle> particles;
  ParticlePars *particle_arr_;
  GPUTrack *gpu_tracks;
  bool external_tracks;
  bool external_geometry;
  bool gpu_tracks_alive;

  void GenerateParticleArray();
  void GenerateParticleIndices(int start, int end);
  void DeleteTracks();

};

} // End namespace na63

#endif /* NA63_SIMULATION_SIMULATOR_H */