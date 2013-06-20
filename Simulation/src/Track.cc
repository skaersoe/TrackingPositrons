#include <cassert>

#include "Simulation/Track.hh"
#include "Simulation/Simulator.hh"
#include "Geometry/Volume.hh"
#include "Geometry/Constants.hh"

namespace na63 {

Float Track::beta() const {
  Float b = momentum_magnitude() / energy();
  assert(b <= 1);
  assert(b == b);
  if (b < 0 || b != b) {
    // Particle is not moving
    b = 0;
  }
  if (b == 1) {
    // Float precision exceeded
    b = 1 - 1e-6;
  }
  assert(b > 0 && b < 1);
  return b;
}

Float Track::gamma() const {
  Float g = Gamma(beta());
  assert(g >= 1);
  return g;
}

ThreeVector Track::beta_vector() const {
  ThreeVector r;
  r = momentum;
  r.Normalize();
  r.Extend(beta());
  return r;
}

void Track::Stop() {
  alive = false;
}

void Track::Kill() {
  momentum[3] = mass();
  Stop();
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
    particle->Query(this,volume->material,dl);
  }
  // Suicide if energy goes below mass
  if (energy() < mass()) Kill();
}

void Track::SpawnChild(Track child) {
  child.mother = this;
  child.volume = volume;
  child.initial_index = initial_index;
  simulator->AddTrackLive(child);
}

void Track::SpawnChildren(std::vector<Track> children) {
  bool live_one = false;
  for (int i=0;i<children.size();i++) {
    SpawnChild(children[i]);
  }
}

Track& Track::operator=(const GPUTrack& gpu_track) {
  alive = (gpu_track.state == ALIVE) ? true : false;
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
  assert(change_energy == change_energy);
  if (change_energy < 0 && simulator->record_energyloss) {
    simulator->FillEnergyLoss(position,-change_energy);
  }
  momentum[3] += change_energy;
  if (momentum[3] < mass()) Kill();
  Float momentum_new = momentum_magnitude();
  momentum.Normalize();
  momentum.Extend(momentum_new);
}

void Track::UpdateMomentum(const FourVector change) {
  momentum -= change;
  if (momentum[3] < mass()) Kill();
}

void Track::SetMomentum(const FourVector momentum_new) {
  momentum = momentum_new;
  if (momentum[3] < mass()) Kill();
}

} // End namespace na63