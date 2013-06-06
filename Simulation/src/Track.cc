#include <cassert>

#include "Simulation/Track.hh"
#include "Simulation/Simulator.hh"
#include "Geometry/Volume.hh"
#include "Geometry/Constants.hh"

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
  Float bp = bx*momentum[0] + by*momentum[1] + bz*momentum[2];
  Float gamma2 = b2 > 0 ? (gamma - 1.0) / b2 : 0.0;

  momentum[0] += gamma2*bp*bx + gamma*bx*momentum[3];
  momentum[1] += gamma2*bp*by + gamma*by*momentum[3];
  momentum[2] += gamma2*bp*bz + gamma*bz*momentum[3];
  momentum[3] = gamma * (momentum[3] + bp);

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
  alive = (gpu_track.state == STATE_ALIVE) ? true : false;
  particle_id = gpu_track.particle_id;
  initial_index = gpu_track.initial_index;
  position = gpu_track.position;
  momentum = gpu_track.momentum;
  vertex_ = gpu_track.vertex;
  return *this;
}

std::ostream& operator<<(std::ostream& os, const Track& t) {
  os << t.particle_id
     << ", " << t.position << ", "
     << t.momentum;
  return os;
}

void Track::UpdateEnergy(const Float change_energy) {
  // p_i' = p_new / p_old * p_i
  momentum[3] += change_energy;
  if (momentum[3] < mass()) Kill();
  Float momentum_new = momentum_magnitude();
  momentum.Normalize();
  momentum.Extend(momentum_new);
}

} // End namespace na63