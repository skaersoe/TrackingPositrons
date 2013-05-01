#include <iostream>

#include "Simulation/Simulator.hh"
#include "Simulation/PropagateGPU.cuh"
#include "Geometry/Library.hh"

namespace na63 {

  Simulator::Simulator(void)
      : Simulator::Simulator(nullptr) {
    geometry = new Geometry();
    external_geometry = false;
  }
  Simulator::Simulator(Geometry *g) {
    geometry = g;
    device = GPU;
    debug = false;
    n_tracks_ = 0;
    external_tracks = true;
    external_geometry = true;
    step_size_ = 0.1;
  }
  Simulator::~Simulator() {
    if (!external_geometry) delete geometry;
    DeleteTracks();
  }
  
  void Simulator::DeleteTracks() {
    if (!external_tracks) {
      delete tracks;
      n_tracks_ = 0;
    }
  }

  int Simulator::SetTracks(Track *t, const unsigned N) {
    DeleteTracks();
    n_tracks_ = N;
    tracks = t;
    for (int i=0;i<N;i++) {
      tracks[i].particle = geometry->GetParticle(tracks[i].particle_id);
      if (tracks[i].particle == nullptr) {
        std::cerr << "Particle with id " << tracks[i].particle_id
                  << " not found." << std::endl;
        return -1;
      }
    }
    return 0;
  }

  Track Simulator::GetTrack(unsigned index) {
    return tracks[index];
  }

  GPUTrack* Simulator::GPUTracks() {
    if (gpu_tracks != nullptr) delete gpu_tracks;
    gpu_tracks = new GPUTrack[n_tracks_];
    for (int i=0;i<n_tracks_;i++) {
      gpu_tracks[i] = tracks[i].GPU();
    }
    GenerateParticleIndices(0,n_tracks_);
    return gpu_tracks;
  }

  void Simulator::CopyBackTracks() {
    for (int i=0;i<n_tracks_;i++) {
      tracks[i] = gpu_tracks[i];
    }
    delete gpu_tracks;
  }

  int Simulator::GenerateParticleIndices(int start, int end) {
    if (gpu_tracks == nullptr) return -1;
    int previous_id = -1;
    int previous_index = -1;
    for (int i=start;i<end;i++) {
      if (gpu_tracks[i].particle_id == previous_id) {
        gpu_tracks[i].particle_index = previous_index;
      } else {
        previous_id = gpu_tracks[i].particle_id;
        int retval = geometry->GetParticleIndex(gpu_tracks[i].particle_id);
        if (retval < 0) {
          std::cerr << "Particle with id " << gpu_tracks[i].particle_id << " not found."
                    << std::endl;
          return -1;
        }
        gpu_tracks[i].particle_index = retval;
        previous_index = gpu_tracks[i].particle_index; 
      }
    }
  }

  void Simulator::Propagate() {

    // Propagate on GPU
    if (device == GPU) {
      geometry->GenerateParameterArrays();
      na63::PropagateGPU(this);
      return;
    }

    // Propagate on CPU
    for (int i=0;i<n_tracks_;i++) {
      Track *t = &tracks[i];
      while (geometry->InBounds(t)) {
        t->Step(step_size_);
        geometry->Query(t);
      }
    }

  }

}