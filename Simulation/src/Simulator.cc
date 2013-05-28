#include <iostream>
#include <thread>
#include <cassert>

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
  external_geometry = true;
  step_size = 0.1;
  particle_arr_ = nullptr;
  cpu_threads = 1;
}
Simulator::~Simulator() {
  if (!external_geometry) delete geometry;
  if (particle_arr_ != nullptr) delete particle_arr_;
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

void Simulator::AddTrack(Track t) {
  t.particle = GetParticle(t.particle_id);
  assert(t.particle != nullptr);
  t.simulator = this;
  tracks.push_back(t);
}

void Simulator::AddTracks(std::vector<Track> t) {
  int previous_id = -1;
  Particle* previous_particle = nullptr;
  for (int i=0;i<t.size();i++) {
    t[i].simulator = this;
    t[i].initial_index = tracks.size() + i;
    int current_id = t[i].particle_id;
    if (current_id == previous_id && current_id != -1) {
      t[i].particle = previous_particle;
    } else {
      previous_particle = t[i].particle = GetParticle(current_id);
      assert(t[i].particle != nullptr);
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

std::vector<Track> Simulator::GetTracks() const {
  return tracks;
}

void Simulator::ClearTracks() {
  tracks.clear();
  std::cout << "All tracks cleared." << std::endl;
}

void Simulator::GPUTracks(GPUTrack* dst) {
  for (int i=0;i<tracks.size();i++) {
    dst[i] = tracks[i].GPU();
    dst[i].initial_index = i;
  }
}

void Simulator::CopyBackTracks(GPUTrack* src, int N) {
  tracks.clear();
  for (int i=0;i<N;i++) {
    Track t;
    t = src[i];
    tracks.push_back(t);
  }
}

void PropagateTrack(Geometry &geometry, Track &track, const Float dl) {
  bool hit = false;
  while (track.alive && geometry.InBounds(track)) {
    track.Step(dl);
    geometry.Query(track);
    if (hit == false && track.volume != nullptr) {
      hit = true;
      // std::cout << "A particle hit the box." << std::endl;
    }
  }
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
  if (cpu_threads > 1) {
    if (debug) std::cout << "Running simulation on " << cpu_threads
        << " threads." << std::endl;
    int i=0;
    int progress = tracks.size()/10;
    std::thread threads[cpu_threads];
    while (i < tracks.size()) {
      int j=0;
      while (j < cpu_threads) {
        // std::cout << "Launching new thread" << std::endl;
        threads[j++] = std::thread(na63::PropagateTrack,std::ref(*geometry),
            std::ref(tracks[i]),step_size);
        if (++i >= tracks.size()) break;
      }
      for (int k=0;k<j;k++) {
        threads[k].join();
      }
      if (debug && i > progress) {
        std::cout << "Propagated " << i << "/" << tracks.size() << " particles..." << std::endl;
        progress += tracks.size()/10;
      }
    }
  } else {
    for (int i=0;i<tracks.size();i++) {
      na63::PropagateTrack(*geometry,tracks[i],step_size);
    }
  }

}

} // End namespace na63