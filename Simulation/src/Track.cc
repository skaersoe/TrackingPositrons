#include "Simulation/Track.hh"
#include "Geometry/Volume.hh"

namespace na63 {
	
void Track::Step(const Float dl) {
	if (!alive_) return;
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
  if (energy() < mass()) Kill();
}

Track& Track::operator=(const GPUTrack& gpu_track) {
  particle_id = gpu_track.particle_id;
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