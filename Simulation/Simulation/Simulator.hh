#ifndef NA63_SIMULATION_SIMULATOR_H
#define NA63_SIMULATION_SIMULATOR_H

#include <thrust/host_vector.h>

#include "Simulation/Track.hh"
#include "Simulation/Particle.hh"
#include "Geometry/Geometry.hh"

namespace na63 {

  typedef enum {
    CPU,
    GPU
  } SimulatorDevice;

  typedef struct {
    SimulatorDevice device; // Run on CPU or GPU
    bool debug;             // Display debugging messages
    bool render;            // Render visually (NYI)
    unsigned N;             // Number of tracks
    Geometry *geometry;
  } SimulatorPars;

  class Simulator {

    public:
      Simulator(void);
      Simulator(SimulatorDevice device, bool debug, unsigned N);

      ~Simulator();

      SimulatorPars pars() {
        return pars_;
      }

      void set_pars(SimulatorPars p) {
        pars_ = p;
      }
      void SetTracks(Track *t, const unsigned N);

      void GenerateTracks();
      void Propagate();

    private:
      SimulatorPars pars_;
      Track *tracks_;
      bool external_tracks;
      Geometry geometry;

      void DeleteTracks();

  };

}

#endif