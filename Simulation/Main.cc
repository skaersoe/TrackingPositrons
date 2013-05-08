#include <iostream>
#include <time.h>
#include <cassert>

#include "Simulation/Simulator.hh"
#include "Geometry/Geometry.hh"
#include "Geometry/Box.hh"
#include "Geometry/Sphere.hh"

#include "Simulation/BetheEnergyLoss.hh"

using namespace na63;

int main(int argc,char *argv[]) {

  // Create geometry
  Geometry geometry;
  geometry.AddMaterial(Material("vacuum",0,0));
  geometry.AddMaterial(Material("iron",26,286));
  geometry.SetBounds(Box("vacuum",ThreeVector(2e2,0,0),ThreeVector(4e2,4e2,4e2)));
  geometry.AddVolume(Box("iron",ThreeVector(2e2,0,0),ThreeVector(1e2,1e2,1e2)));
  //geometry.SetBounds(Box("vacuum",ThreeVector(2e2,0,0),ThreeVector(4e2,4e2,4e2)));
  //geometry.AddVolume(Box("iron",ThreeVector(2e2,0,0),ThreeVector(1e2,1e2,1e2)));

  // Create Simulator object
  Simulator sim = Simulator(&geometry);
  Particle muon = Particle("muon",13,-1,105.6583715);
  muon.RegisterProcess(BetheEnergyLoss);
  sim.AddParticle(Particle("electron",11,-1,0.510998910));
  sim.AddParticle(muon);

  // Set some default values
  sim.device = GPU;
  sim.debug = false;
  unsigned N = 0;

  // Parse input
  if (argc < 2) {
    std::cerr << "Missing input: number of particles to simulate." << std::endl;
    return -1;
  }
  sscanf(argv[1],"%d",&N);
  if (N <= 0) {
    std::cerr << "Invalid number of particles." << std::endl;
    return -1;
  }
  for (int i=2;i<argc;i++) {
    std::string token(argv[i]);
    if (token == "-CPU")
      sim.device = CPU;
    else if (token == "-debug")
      sim.debug = true;
    else {
      std::cout << "Unrecognized argument: " << token << std::endl;
      return -1;
    }
  } // Finished parsing input arguments

  std::cout << "Simulator starting." << std::endl;

  // Tracks
  std::cout << "Generating tracks..." << std::endl;
  std::vector<Track> t;
  const Float arc = kPi/4;
  const Float v = 0.95; // % of c
  const Float m = 105.6583715;
  const Float E = Gamma(v) * m;
  for (int i=0;i<N;i++) {
    Float angle = -arc + 2*arc * ((Float)i / (Float)N);
    Float vx = v * cos(angle);
    Float vy = v * sin(angle);
    Float px = Gamma(vx) * m * vx;
    Float py = Gamma(vy) * m * vy;
    Float pz = 0;
    t.push_back(Track(13,FourVector(),FourVector(px,py,pz,E)));
  }
  sim.AddTracks(t);

  // Propagate
  clock_t timer = clock();
  sim.Propagate();
  timer = clock() - timer;
  std::cout << "Propagated in " << (float)timer/CLOCKS_PER_SEC << " seconds." << std::endl;

  std::cout << "Simulator exiting." << std::endl;
  return 0;
}