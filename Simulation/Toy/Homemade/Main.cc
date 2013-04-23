#include <iostream>
#include <cmath>
#include <vector>

#include "Track.hh"
#include "SimpleBox.hh" 
#include "Geometry.hh"

using namespace std;

int hits = 0;

float Gamma(float v) {
  return 1/sqrt(1-pow(v,2));
}

void Propagate(Track *track, Geometry *geometry) {
  while (geometry->InBounds(track)) {
    track->step(0.01);
    if (geometry->BoundaryDetection(track) != NULL) {
      hits++;
      return;
    }
  }
}

void PrintTrack(int i, Track track) {
  cout << "Track " << i << ": (" << track.x() << "," << track.y() << "," << track.z() << ")" << endl;
}

int main(void) {

  // Parameters
  const int particle_count = 256;
  const float arc = M_PI/4;
  const float m = 0.510998910; // MeV/c^2
  const float v = 0.95; // % of c
  const float E = Gamma(v) * m;

  // Tracks
  cout << "Generating particles..." << endl;
  vector<Track> particles;
  for (int i=0;i<particle_count;i++) {
    float angle = -arc + 2*arc * ((float)i / (float)particle_count);
    float vx = v * cos(angle);
    float vy = v * sin(angle);
    float px = Gamma(vx) * m * vx;
    float py = Gamma(vy) * m * vy;
    float pz = 0;
    particles.push_back(Track(ThreeVector(0,0,0),FourMomentum(ThreeVector(px,py,pz),E),-1));
  }

  // Experiment
  Geometry geometry(SimpleBox(ThreeVector(2e3,0,0),ThreeVector(4e3,4e3,0)));
  geometry.AddBox(SimpleBox(ThreeVector(2e3,0,0),ThreeVector(1e3,1e3,0)));

  // Run
  cout << "Propagating particles..." << endl;
  for (int i=0;i<particle_count;i++)
    Propagate(&particles[i],&geometry);

  cout << hits << "/" << particle_count << " particles hit the iron box." << endl;

  return 0;
}