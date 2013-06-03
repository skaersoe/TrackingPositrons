#include <iostream>
#include <vector>
#include "Geometry/Constants.hh"
#include "Simulation/Simulator.hh"
#include "Simulation/Particle.hh"
#include "Simulation/Track.hh"
#include "Geometry/Geometry.hh"
#include "Geometry/Sphere.hh"
#include <TRandom3.h>

using namespace na63;

int main(void) {
  
  const int n_tracks = 32;

  Geometry geometry;
  geometry.AddMaterial(Material("vacuum",0.0,0.0,0.0));
  geometry.AddMaterial(Material("iron",26.0,286.0,13.84));
  geometry.SetBounds(Sphere("vacuum",ThreeVector(0,0,0),3e2));
  geometry.AddVolume(Sphere("iron",ThreeVector(0,0,0),1e2));

  Simulator simulator(&geometry);
  simulator.device = GPU;
  simulator.debug = true;
  simulator.step_size = 0.1;
  simulator.sorting = RADIUS;
  simulator.pool_size = n_tracks*20;
  simulator.AddParticle(Particle("photon",22,0));

  const Float E = 10 * MeV;

  std::vector<Track> tracks;
  for (int i=0;i<n_tracks;i++) {
    Float theta = 2 * kPi * gRandom->Uniform();
    Float phi = kPi * gRandom->Uniform();
    ThreeVector direction = SphericalToCartesian(E,phi,theta);
    tracks.push_back(Track(22,0,FourVector(),FourVector(direction,E)));
  }
  simulator.AddTracks(tracks);

  simulator.Propagate(); 
  tracks = simulator.GetTracks();

  for (int i=0;i<tracks.size();i++) {
    if (tracks[i].alive) {
      std::cout << "Track " << i << " did not finish." << std::endl;
    }
  }

  std::cout << tracks.size() << " tracks were returned." << std::endl;

  return 0;
}