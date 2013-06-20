#include "Simulation/Particle.hh"
#include "Simulation/Track.hh"

namespace na63 {

Particle::Particle(const char* n, int id, Float mass) : name_(n) {
  id_ = id;
  mass_ = mass;
  index = -1;
}

ParticlePars Particle::GPU() const {
  ParticlePars retval;
  retval.id = id_;
  retval.mass = mass_;
  return retval;
}

void Particle::RegisterProcess(Process *process) {
  processes.push_back(process);
}

void Particle::Query(Track* track, const Material* material, const Float dl) const {
  for (int i=0;i<processes.size();i++) {
    if (!track->alive) return;
    processes[i]->Query(track,material,dl);
  }
}

} // End namespace na63