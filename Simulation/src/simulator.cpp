#include <iostream>

#include "Simulation/Simulator.h"
#include "Simulation/GenerateTracks.h"
#include "Simulation/Propagate.h"

namespace na63 {

  Simulator::Simulator(void)
    : Simulator::Simulator(GPU, false) {}
  Simulator::Simulator(SimulatorDevice device, bool debug) {
    pars_ = {
      .device = device,
      .debug = debug,
      .render = false,
      .N = 0,
      .geometry = &geometry
    };
    external_tracks = true;
  }
  Simulator::~Simulator() {
    DeleteTracks();
  }
  
  void Simulator::DeleteTracks() {
    if (!external_tracks) {
      delete tracks_;
      pars_.N = 0;
    }
  }

  void Simulator::SetTracks(Track *t, const unsigned N) {
    DeleteTracks();
    pars_.N = N;
    tracks_ = t;
  }

  /**
   * Generates some hardcoded electrons.
   */
  void Simulator::GenerateTracks() {
    DeleteTracks();
    tracks_ = new Track[pars_.N];
    na63::GenerateTracks(tracks_,pars_);
    external_tracks = false;
  }

  void Simulator::Propagate() {
    geometry.GenerateParameterArrays();
    na63::Propagate(tracks_,pars_);
  }

}