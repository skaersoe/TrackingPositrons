#ifndef NA63_SIMULATION_SIMULATOR_H
#define NA63_SIMULATION_SIMULATOR_H

#include <thrust/host_vector.h>

#include "Simulation/Track.h"
#include "Simulation/Arguments.h"

namespace na63 {

  class Simulator {

    public:
      Simulator(void);
      Simulator(SimulatorPars pars);

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

      void DeleteTracks();

  };

}

#endif