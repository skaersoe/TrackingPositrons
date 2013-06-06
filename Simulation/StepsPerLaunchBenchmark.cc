#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>

// ROOT crap
#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>
#include "TRandom3.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"

#include "Simulation/Simulator.hh"
#include "Geometry/Geometry.hh"
#include "Geometry/Box.hh"
#include "Geometry/Sphere.hh"
#include "Simulation/GetTime.hh"

#include "Simulation/BetheEnergyLoss.hh"

using namespace na63;

int main(void) {

  const int number_of_runs = 10;
  const int number_of_tracks = 4096;
  int step_runs[] = {
    50,
    100,
    150,
    200,
    400,
    600,
    800,
    1000,
    1500,
    2000
  };
  double *benchmarks = new double[number_of_runs];

  TRandom3 rng;

  // Set up geometry
  Geometry geometry;
  geometry.AddMaterial(Material("vacuum",0,0,0,0));
  geometry.AddMaterial(Material("iron",kIronAtomicNumber,kIronDensity,kIronMeanExcitationPotential,kIronRadiationLength));
  geometry.SetBounds(Sphere("vacuum",ThreeVector(0,0,0),1e10));
  geometry.AddVolume(Sphere("iron",ThreeVector(0,0,0),1e10));

  // Create simulator
  Simulator simulator = Simulator(&geometry);
  Particle muon = Particle("muon",13.0,kMuonMass);
  muon.RegisterProcess(BetheEnergyLoss);
  simulator.AddParticle(muon);
  simulator.step_size = 0.1;

  // Set parameters
  simulator.device = GPU;
  simulator.debug = false;

  const Float energy = 1e5; // 100 GeV
  std::vector<Track> tracks;
  for (int i=0;i<number_of_tracks;i++) {
    tracks.push_back(Track(13,-1,FourVector(),
            FourVector(energy/kMuonMass,0,0,energy)));
  }

  for (int i=0;i<number_of_runs;i++) {

    simulator.steps_per_launch = step_runs[i];

    // Propagate on GPU
    simulator.AddTracks(tracks);
    double elapsed = InSeconds(GetTime());
    std::cout << "Running with step size " << step_runs[i] << "..." << std::endl;
    simulator.Propagate();
    benchmarks[i] = InSeconds(GetTime()) - elapsed;
    simulator.ClearTracks();

  }

  std::ofstream outfile;
  outfile.open("Data/stepsperlaunch_benchmark");
  outfile << "steps,time in seconds" << std::endl;
  for (int i=0;i<number_of_runs;i++) {
    outfile << step_runs[i] << "," << benchmarks[i] << std::endl;
  }
  outfile.close();

  return 0;
}