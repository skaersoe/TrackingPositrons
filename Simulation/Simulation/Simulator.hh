#ifndef NA63_SIMULATION_SIMULATOR_H
#define NA63_SIMULATION_SIMULATOR_H

#include <vector>

#include "Simulation/Track.hh"
#include "Simulation/Particle.hh"
#include "Geometry/Geometry.hh"

namespace na63 {

enum Device {CPU,GPU};
enum SortMethod {X,Y,Z,RADIUS};

class Simulator {

public:
  // Allow these parameters to be altered freely
  Device device;
  SortMethod sorting;
  bool debug;
  int cpu_threads;
  int pool_size;
  Float step_size;
  int steps_per_launch;
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
  std::vector<Track> GetTracks() const;
  void AddTrack(Track t);
  void AddTracks(std::vector<Track> t);
  void ClearTracks();
  //void SetTracks(std::vector<Track> t);
  unsigned TrackSize() { return tracks.size(); }
  void GPUTracks(GPUTrack* dst);
  void CopyBackTracks(GPUTrack* src, int N);

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
  bool external_geometry;

  void GenerateParticleArray();
  void GenerateParticleIndices(int start, int end);
  void DeleteTracks();

};

} // End namespace na63

#endif /* NA63_SIMULATION_SIMULATOR_H */