#include <iostream>
#include <thread>
#include <cassert>

#include "Simulation/Simulator.hh"
#include "Simulation/PropagateGPU.cuh"
#include "Geometry/Library.hh"
#include "Simulation/GetTime.hh"

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
  pool_size = 2;
  sorting = RADIUS;
  steps_per_launch = 100;
  thread_multiplier = 1;
  record_energyloss = false;
  gpu_bremsstrahlung = true;
  secondary_threshold = 2.0*kElectronMass;
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
  AddTracksGeneric(t,&tracks);
}

void Simulator::AddTrackLive(Track t) {
  AddTracksGeneric(t,&tracks_waiting);
}

void Simulator::AddTracksLive(std::vector<Track> t) {
  for (int i=0;i<t.size();i++) AddTrackLive(t[i]);
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
  tracks.insert(tracks.end(),t.begin(),t.end());
  std::cout << "Tracks added. New size is " << tracks.size() 
            << "." << std::endl;
}

Track Simulator::GetTrack(int index) const {
  return tracks[index];
}

Track* Simulator::GetTrackAddress(int index) {
  return &tracks[index];
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
  std::cout << "Copying back tracks... ";
  tracks.clear();
  int j=0;
  for (int i=0;i<N;i++) {
    assert(src[i].state != WAITING);
    if (src[i].state == FREE) continue;
    Track t;
    t = src[i];
    tracks.push_back(t);
    j++;
  }
  std::cout << j << " tracks were successfully copied back." << std::endl;
}

void Simulator::PrintTracks() const {
  for (int i=0;i<tracks.size();i++) {
    std::cout << tracks[i].particle_id << ", " << tracks[i].position << ", "
              << tracks[i].momentum << ((!tracks[i].alive) ? " (dead)" : "") << std::endl;
  }
  for (int i=0;i<tracks_waiting.size();i++) {
    std::cout << tracks_waiting[i].particle_id << ", " << tracks_waiting[i].position << ", "
              << tracks_waiting[i].momentum << ((!tracks_waiting[i].alive) ? " (dead)" : "") << std::endl;
  }
}

void Simulator::SetBenchmark(double t) {
  benchmark = t;
}

double Simulator::GetBenchmark() const {
  return benchmark;
}

void Simulator::AppendWaitingTracks() {
  tracks.insert(tracks.end(),tracks_waiting.begin(),tracks_waiting.end());
  tracks_waiting.clear();
}

void Simulator::RecordEnergyLoss(EnergyLoss *energyloss_) {
  record_energyloss = true;
  energyloss = energyloss_;
}

void Simulator::FillEnergyLoss(const Float x, const Float y, const Float z,
    const Float E) {
  assert(E >= 0);
  energyloss->Fill(x,y,z,E);
}

void Simulator::FillEnergyLoss(const FourVector& position, const Float E) {
  FillEnergyLoss(position[0],position[1],position[2],E);
}

void PropagateTrack(Geometry &geometry, Track &track, const Float dl) {
  Float travelled = 0;
  while (track.alive && geometry.InBounds(track)) {
    track.Step(dl);
    geometry.Query(track);
    travelled += dl;
    if (travelled > 1e6) {
      std::cout << "Stopping track that has propagated for too long: " << track << std::endl;
      break;
    }
  }
  track.Stop();
}

void Simulator::Propagate() {

  // Propagate on GPU
  if (device == GPU) {
    std::cout << "Running simulation GPU." << std::endl;
    if (debug) std::cout << "Generating parameter arrays...";
    geometry->GenerateParameterArrays();
    GenerateParticleArray();
    if (debug) std::cout << " OK" << std::endl;
    na63::PropagateGPU(this);
    return;
  }

  // Propagate on CPU
  std::cout << "Running simulation on CPU." << std::endl;
  benchmark = InSeconds(GetTime());
  if (cpu_threads > 1) {
    if (debug) std::cout << "Running simulation on " << cpu_threads
        << " threads." << std::endl;
    int i=0;
    int progress = tracks.size()/10;
    int previous = 0;
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
      AppendWaitingTracks();
      if (debug) {
        if (i > progress) {
          std::cout << "Propagated " << i << "/" << tracks.size() << " particles..." << std::endl;
          previous = i;
          progress += tracks.size()/10;
        } else if (i > previous + 1000) {
          std::cout << "Propagated " << i << "/" << tracks.size() << " particles..." << std::endl;
          previous = i;
          progress += tracks.size()/10;
        }
      }
    }
  } else {
    for (int i=0;i<tracks.size();i++) {
      na63::PropagateTrack(*geometry,tracks[i],step_size);
    }
  }
  benchmark = InSeconds(GetTime()) - benchmark;
  std::cout << "Propagation ran for " << benchmark << " seconds, and returned " << TrackSize() << " tracks." << std::endl;

}

} // End namespace na63