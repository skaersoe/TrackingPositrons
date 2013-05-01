#include "Track.hh"
#include "Physics.hh"

void Track::Step(Float dl) {
  ThreeVector change = SphericalToCartesian(dl,momentum.Theta(),momentum.Phi());
  position[0] += change[0];
  position[1] += change[1];
  position[2] += change[2];
  position[3] += dl * (E() / (p() * c));
  BetheIron(this,dl);
  return;
}