#ifndef NA63_SIMULATION_SIMULATOR_H
#define NA63_SIMULATION_SIMULATOR_H

#include <vector>

#include "Simulation/Track.hh"
#include "Simulation/Particle.hh"
#include "Geometry/Geometry.hh"
#ifndef __CUDACC__
#include "Simulation/EnergyLoss.hh"
#else
class EnergyLoss;
#endif
namespace na63 {

enum Device {CPU,GPU};
enum SortMethod {X,Y,Z,RADIUS,PARTICLE};

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
  int thread_multiplier;
  Geometry *geometry;
  bool record_energyloss;
  bool gpu_bremsstrahlung;
  Float secondary_threshold;

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
  Track* GetTrackAddress(int index);
  std::vector<Track> GetTracks() const;
  void AddTrack(Track t);
  void AddTrackLive(Track t);
  void AddTracks(std::vector<Track> t);
  void AddTracksLive(std::vector<Track> t);
  void ClearTracks();
  //void SetTracks(std::vector<Track> t);
  unsigned TrackSize() { return tracks.size() + tracks_waiting.size(); }
  void GPUTracks(GPUTrack* dst);
  void CopyBackTracks(GPUTrack* src, int N);
  void PrintTracks() const;
  void RecordEnergyLoss(EnergyLoss *energyloss_);
  void FillEnergyLoss(const Float x, const Float y, const Float z,
      const Float E);
  void FillEnergyLoss(const FourVector& position, const Float E);

  /**
   * Generates some hardcoded electrons.
   */
  // void GenerateTracks();
  /**
   * Runs the propagation.
   */
  void Propagate();
  void PropagateTrack(Track& track);
  void SetBenchmark(double t);
  double GetBenchmark() const;

private:
  std::vector<Track> tracks;
  std::vector<Track> tracks_waiting;
  std::vector<Particle> particles;
  ParticlePars *particle_arr_;
  bool external_geometry;
  double benchmark;
  EnergyLoss *energyloss;

  void GenerateParticleArray();
  void GenerateParticleIndices(int start, int end);
  void DeleteTracks();
  void AppendWaitingTracks();

  inline void AddTracksGeneric(Track t,std::vector<Track> *tgt) {
    t.particle = GetParticle(t.particle_id);
    assert(t.particle != nullptr);
    if (t.energy() <= t.mass()) t.Kill();
    t.simulator = this;
    tgt->push_back(t);
  }

};

} // End namespace na63

#endif /* NA63_SIMULATION_SIMULATOR_H */