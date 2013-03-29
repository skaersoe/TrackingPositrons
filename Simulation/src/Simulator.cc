#include <iostream>

#include "Simulation/Simulator.hh"
#include "Simulation/GenerateTracks.cuh"
#include "Simulation/Propagate.cuh"

namespace na63 {

  Simulator::Simulator(void)
    : Simulator::Simulator(GPU, false, 0) {}
  Simulator::Simulator(SimulatorDevice device, bool debug, unsigned N) {
    pars_ = {
      .device = device,
      .debug = debug,
      .render = false,
      .N = N,
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