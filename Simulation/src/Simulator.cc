#include <iostream>

#include "Simulation/Simulator.hh"
#include "Simulation/GenerateTracks.cuh"
#include "Simulation/Propagate.cuh"

namespace na63 {

  Simulator::Simulator(void)
      : Simulator::Simulator(nullptr) {
    geometry = new Geometry();
    pars.geometry = geometry;
    external_geometry = false;
  }
  Simulator::Simulator(Geometry *g) {
    geometry = g;
    pars = {
      .device = GPU,
      .debug = false,
      .render = false,
      .N = 0,
      .geometry = geometry
    };
    external_tracks = true;
    external_geometry = true;
  }
  Simulator::~Simulator() {
    if (!external_geometry) delete geometry;
    DeleteTracks();
  }
  
  void Simulator::DeleteTracks() {
    if (!external_tracks) {
      delete tracks_;
      pars.N = 0;
    }
  }

  void Simulator::SetTracks(Track *t, const unsigned N) {
    DeleteTracks();
    pars.N = N;
    tracks_ = t;
  }

  /**
   * Generates some hardcoded electrons.
   */
  void Simulator::GenerateTracks() {
    DeleteTracks();
    tracks_ = new Track[pars.N];
    na63::GenerateTracks(tracks_,pars);
    external_tracks = false;
  }

  int Simulator::GenerateParticleIndices(int start, int end) {
    int previous_id = -1;
    int previous_index = -1;
    for (int i=start;i<end;i++) {
      if (tracks_[i].particle_id == previous_id) {
        tracks_[i].particle_id = previous_index;
      } else {
        previous_id = tracks_[i].particle_id;
        int retval = geometry->GetParticleIndex(tracks_[i].particle_id);
        if (retval < 0) {
          std::cerr << "Particle with id " << tracks_[i].particle_id << " not found."
                    << std::endl;
          return -1;
        }
        tracks_[i].particle_id = retval;
        previous_index = tracks_[i].particle_id; 
      }
    }
  }

  void Simulator::Propagate() {
    if (GenerateParticleIndices(0,pars.N) < 0) return;
    geometry->GenerateParameterArrays();
    na63::Propagate(tracks_,pars);
  }

}