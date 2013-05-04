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
    external_tracks = true;
    external_geometry = true;
    step_size_ = 0.1;
  }
  Simulator::~Simulator() {
    if (!external_geometry) delete geometry;
  }

  /*void Simulator::SetTracks(std::vector<Track> t) {
    tracks = t;
  }*/

  void Simulator::AddTracks(Track t, const int N) {
    if ((t.particle = geometry->GetParticle(t.particle_id)) == nullptr) {
      std::cerr << "Particle with id " << t.particle_id
                << " not found. Tracks were not added." << std::endl;
    }
    for (int i=0;i<N;i++) {
      tracks.push_back(t);
    }
    std::cout << "Tracks added. New size is " << tracks.size() 
              << "." << std::endl;
  }

  Track Simulator::GetTrack(int index) const {
    return tracks[index];
  }

  GPUTrack* Simulator::GPUTracks() {
    if (gpu_tracks != nullptr) delete gpu_tracks;
    gpu_tracks = new GPUTrack[tracks.size()];
    for (int i=0;i<tracks.size();i++) {
      gpu_tracks[i] = tracks[i].GPU();
    }
    GenerateParticleIndices(0,tracks.size());
    return gpu_tracks;
  }

  void Simulator::CopyBackTracks() {
    for (int i=0;i<tracks.size();i++) {
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
      if (debug) std::cout << "Generating parameter arrays...";
      geometry->GenerateParameterArrays();
      if (debug) std::cout << " OK" << std::endl;
      na63::PropagateGPU(this);
      return;
    }

    // Propagate on CPU
    for (int i=0;i<tracks.size();i++) {
      std::cout << tracks[i] << std::endl;
      while (geometry->InBounds(tracks[i])) {
        tracks[i].Step(step_size_);
        geometry->Query(tracks[i]);
        if (tracks[i].volume != nullptr) {
          std::cout << "Particle " << i << " hit something." << std::endl;
          continue;
        }
      }
    }

  }

}