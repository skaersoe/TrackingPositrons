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
    Geometry *geometry;     // Pointer to associated geometry
  } SimulatorPars;

  class Simulator {

    public:
      // Allow parameters to be altered freely
      SimulatorPars pars;

      Simulator(void);
      Simulator(SimulatorDevice device, bool debug, unsigned N);
      ~Simulator();

      void set_pars(SimulatorPars p) {
        pars = p;
      }
      void SetTracks(Track *t, const unsigned N);

      /**
       * Generates some hardcoded electrons.
       */
      void GenerateTracks();
      /**
       * Runs the propagation.
       */
      void Propagate();

    private:
      Track *tracks_;
      bool external_tracks;
      Geometry geometry;

      void DeleteTracks();

  };

}

#endif /* NA63_SIMULATION_SIMULATOR_H */