#include <iostream>
#include <cmath>
#include <vector>
#include <random>

#include "Track.hh"
#include "SimpleBox.hh" 
#include "Geometry.hh"
#include "Library.hh"

using namespace std;

int hits = 0;

Float Gamma(Float v) {
  return 1/sqrt(1-pow(v,2));
}

void Propagate(Track *track, Geometry *geometry) {
  bool hit = false;
  while (geometry->InBounds(track)) {
    track->Step(0.1);
  }
}

void PrintTrack(int i, Track track) {
  cout << "Track " << i << ": (" << track.x() << "," << track.y() << "," << track.z() << "," << track.t() << ")" << endl;
}

int main(void) {

  // Parameters
  const int particle_count = 1;
  const Float arc = M_PI/4;
  const Float m = 0.510998910; // MeV/c^2
  const Float v = 0.95; // % of c
  const Float E = Gamma(v) * m;

  // Tracks
  cout << "Generating particles..." << endl;
  vector<Track> particles;
  for (int i=0;i<particle_count;i++) {
    Float angle = -arc + 2*arc * ((Float)i / (Float)particle_count);
    Float vx = v * cos(angle);
    Float vy = v * sin(angle);
    Float px = Gamma(vx) * m * vx;
    Float py = Gamma(vy) * m * vy;
    Float pz = 0;
    particles.push_back(Track(FourVector(),FourVector(px,py,pz,E),-1,m));
  }

  // Experiment
  Geometry geometry(SimpleBox(ThreeVector(2e3,0,0),ThreeVector(4e3,4e3,4e3)));
  geometry.AddBox(SimpleBox(ThreeVector(2e3,0,0),ThreeVector(1e3,1e3,1e3)));

  // Run
  cout << "Propagating particles..." << endl;
  for (int i=0;i<particle_count;i++)
    Propagate(&particles[i],&geometry);

  cout << hits << "/" << particle_count << " particles hit the iron box." << endl;

  return 0;
}