#include <iostream>

#include "Simulation/Simulator.h"
#include "Simulation/GenerateTracks.h"
#include "Simulation/Propagate.h"

namespace na63 {

  Simulator::Simulator(void)
    : Simulator::Simulator({ .device = GPU, .debug = false }) {}
  Simulator::Simulator(SimulatorPars p) {
    set_pars(p);
    tracks_ = NULL;
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

  void Simulator::GenerateTracks() {
    DeleteTracks();
    tracks_ = new Track[pars_.N];
    na63::GenerateTracks(tracks_,pars_);
    external_tracks = false;
  }

  void Simulator::Propagate() {
    na63::Propagate(tracks_,pars_);
  }

}