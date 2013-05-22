#include "Simulation/Track.hh"
#include "Geometry/Volume.hh"
#include "Simulation/Simulator.hh"
#include <cassert>

namespace na63 {

ThreeVector Track::beta_vector() const {
  ThreeVector r;
  r = momentum;
  r.Normalize();
  r.Extend(beta());
  return r;
}

void Track::Boost(const Float bx, const Float by, const Float bz) {

  Float b2 = pow(bx,2) + pow(by,2) + pow(bz,2);

  Float gamma = 1.0 / sqrt(1.0 - b2);
  Float bp = bx*position[0] + by*position[1] + bz*position[2];
  Float gamma2 = b2 > 0 ? (gamma - 1.0) / b2 : 0.0;

  position[0] += gamma2*bp*bx + gamma*bx*position[3];
  position[1] += gamma2*bp*by + gamma*by*position[3];
  position[2] += gamma2*bp*bz + gamma*bz*position[3];
  position[3] = gamma * (position[3] + bp);

}
	
void Track::Step(const Float dl) {
	if (!alive) return;
	ThreeVector change = SphericalToCartesian(
	    dl,momentum.Theta(),momentum.Phi());
	position[0] += change[0];
	position[1] += change[1];
	position[2] += change[2];
	position[3] += dl * (energy() / (momentum_magnitude() * kC));
	// Let's do physics!
	if (volume != nullptr) {
    particle->Query(*this,*volume->material,dl);
  }
  // Suicide if energy goes below mass
  if (energy() < mass()) Kill();
}

void Track::SpawnChild(Track child) {
  child.mother = this;
  child.volume = volume;
  child.initial_index = initial_index;
  simulator->AddTrack(child);
}

Track& Track::operator=(const GPUTrack& gpu_track) {
  particle_id = gpu_track.particle_id;
  initial_index = gpu_track.initial_index;
  position = gpu_track.position;
  momentum = gpu_track.momentum;
  return *this;
}

std::ostream& operator<<(std::ostream& os, const Track& t) {
  if (t.particle != nullptr) {
    os << t.particle->name();
  } else {
    os << t.particle_id;
  }
  os << ", " << t.position << ", "
     << t.momentum;
  return os;
}

} // End namespace na63