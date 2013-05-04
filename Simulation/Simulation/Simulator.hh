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
      Geometry *geometry;

      Simulator(void);
      Simulator(Geometry *geometry);
      ~Simulator();

      Track GetTrack(int index) const;
      void AddTracks(Track t, const int N);
      //void SetTracks(std::vector<Track> t);
      unsigned TrackSize() { return tracks.size(); }
      Float step_size() { return step_size_; }
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

    private:
      std::vector<Track> tracks;
      GPUTrack *gpu_tracks;
      bool external_tracks;
      bool external_geometry;
      Float step_size_;

      int GenerateParticleIndices(int start, int end);
      void DeleteTracks();

  };

}

#endif /* NA63_SIMULATION_SIMULATOR_H */