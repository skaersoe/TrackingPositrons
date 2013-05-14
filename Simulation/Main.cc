#include <iostream>
#include <cassert>

// ROOT crap
#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>
#include "TRandom2.h"
#include "TH1.h"
#include "TH2.h"

#include "Simulation/Simulator.hh"
#include "Geometry/Geometry.hh"
#include "Geometry/Box.hh"
#include "Geometry/Sphere.hh"
#include "Simulation/GetTime.hh"

#include "Simulation/BetheEnergyLoss.hh"
#include "Simulation/MultipleScattering.hh"

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
  sim.AddParticle(muon);

  // Set some default values
  bool plot = false;
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
    else if (token == "-plot")
      plot = true;
    else if (token == "-threads") {
      sscanf(argv[++i],"%d",&sim.cpu_threads);
      if (sim.cpu_threads <= 0) {
        std::cerr << "Invalid number of threads." << std::endl;
        return -1;
      }
    }
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
  const Float kMuonMass = 105.6583715;
  const Float E = Gamma(v) * kMuonMass;
  Float angles[N];
  TRandom3 rng((size_t)clock());
  for (int i=0;i<N;i++) {
    Float angle = arc * (rng.Gaus(0,1));
    angles[i] = angle;
    Float vx = v * cos(angle);
    Float vy = v * sin(angle);
    Float px = Gamma(vx) * kMuonMass * vx;
    Float py = Gamma(vy) * kMuonMass * vy;
    Float pz = 0;
    t.push_back(Track(13,FourVector(),FourVector(px,py,pz,E)));
  }
  sim.AddTracks(t);

  // Propagate
  timespec before, after;
  before = GetTime();
  sim.Propagate();
  after = GetTime();
  double elapsed = after.tv_sec - before.tv_sec;
  std::cout << "Propagated in " << elapsed << " seconds." << std::endl;

  // Do ROOT crap
  if (plot) {
    TApplication app("app",&argc,argv);
    TCanvas canvas;
    /*TH2F hist_final_position(
      "Final position","Final position;x;y",
      100,0,4e2+5,100,-2e2-5,2e2+5
    );*/
    TGraph graph_energy_angle;
    for (int i=0;i<N;i++) {
      Track t = sim.GetTrack(i);
      // std::cout << "Died at " << t.position << std::endl;
      //hist_final_position.Fill(t.position[0],t.position[1]);
      graph_energy_angle.SetPoint(i,angles[i],t.energy());
    }
    //hist_final_position.Draw("COLZ");
    graph_energy_angle.Draw("same ap");
    gPad->WaitPrimitive();
  }

  std::cout << "Simulator exiting." << std::endl;
  return 0;
}