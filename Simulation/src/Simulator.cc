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
  particle_arr_    = nullptr;
}
Simulator::~Simulator() {
  if (!external_geometry) delete geometry;
  delete particle_arr_;
}

int Simulator::GetParticleIndex(int id) const {
  for (int i=0;i<particles.size();i++) {
    if (particles[i].id() == id) {
      return i;
    }
  }
  return -1;
}

Particle* Simulator::GetParticle(int id) {
  int index = GetParticleIndex(id);
  if (index == -1) return nullptr;
  return &particles[GetParticleIndex(id)];
}

ParticlePars* Simulator::particle_arr() {
  if (particle_arr_ == nullptr) {
    GenerateParticleArray();
  }
  return particle_arr_;
}

void Simulator::GenerateParticleArray() {
  ParameterVectorToArray(&particles,&particle_arr_);
}

void Simulator::AddTracks(std::vector<Track> t) {
  int previous_id = -1;
  Particle* previous_particle = nullptr;
  for (int i=0;i<t.size();i++) {
    int current_id = t[i].particle_id;
    if (current_id == previous_id && current_id != -1) {
      t[i].particle = previous_particle;
    } else {
      previous_particle = t[i].particle = GetParticle(current_id);
      previous_id = current_id;
    }
  }
  tracks.assign(t.begin(),t.end());
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
  return gpu_tracks;
}

void Simulator::CopyBackTracks() {
  for (int i=0;i<tracks.size();i++) {
    tracks[i] = gpu_tracks[i];
  }
  delete gpu_tracks;
}

void Simulator::Propagate() {

  // Propagate on GPU
  if (device == GPU) {
    if (debug) std::cout << "Generating parameter arrays...";
    geometry->GenerateParameterArrays();
    GenerateParticleArray();
    if (debug) std::cout << " OK" << std::endl;
    na63::PropagateGPU(this);
    return;
  }

  // Propagate on CPU
  for (int i=0;i<tracks.size();i++) {
    bool hit = false;
    // std::cout << tracks[i] << std::endl;
    while (geometry->InBounds(tracks[i])) {
      tracks[i].Step(step_size_);
      geometry->Query(tracks[i]);
      if (tracks[i].volume != nullptr && hit == false) {
        hit = true;
        std::cout << "Particle " << i << " hit something!" << std::endl;
      }
    }
  }

}

} // End namespace na63