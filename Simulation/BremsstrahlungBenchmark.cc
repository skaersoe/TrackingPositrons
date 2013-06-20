#include <iostream>
#include <fstream>
#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <TRandom3.h>
#include <cassert>

#include "Geometry/Geometry.hh"
#include "Geometry/SimpleBox.hh"
#include "Simulation/GEANT4Bremsstrahlung.hh"
#include "Simulation/GEANT4PairProduction.hh"
#include "Simulation/BetheEnergyLoss.hh"
#include "Simulation/Simulator.hh"
#include "Simulation/GetTime.hh"

using namespace na63;

int main(int argc,char *argv[]) {

  const unsigned n_tracks = 512;
  const unsigned runs = 18;
  const unsigned steps[] = {100,200,300,400,600,800,900,1000,1100,1200,1500,2000,2500,3000,3500,4000,4500,5000};
  double *benchmark = new double[runs + 1];

  // Run on CPU
  bool runcpu = false;
  for (int i=1;i<argc;i++) {
    std::string arg(argv[i]);
    if (arg == "CPU") {
      runcpu = true;
      break;
    }
  }


  // Geometry
  Geometry geometry;
  const Float scale = 1e3;
  geometry.AddMaterial(Material("vacuum",0.0,0.0,0.0,0.0,0.0));
  geometry.AddMaterial(Material("iron",kIronAtomicNumber,kIronDensity,kIronAtomicWeight,kIronMeanExcitationPotential,kIronRadiationLength));
  geometry.SetBounds(SimpleBox("vacuum",ThreeVector(5*scale,0,0),ThreeVector(10*scale,2*scale,2*scale)));
  // geometry.AddVolume(SimpleBox("iron",ThreeVector(50,1e2+5,0),ThreeVector(30,2e2,4e2)));
  // geometry.AddVolume(SimpleBox("iron",ThreeVector(50,-(1e2+5),0),ThreeVector(30,2e2,4e2)));
  // geometry.AddVolume(SimpleBox("iron",ThreeVector(2*scale,0,0),ThreeVector(1*scale,1*scale,1*scale)));
  // geometry.AddVolume(SimpleBox("iron",ThreeVector(1*scale,1*scale,1*scale),ThreeVector(1*scale,1*scale,1*scale)));
  // geometry.AddVolume(SimpleBox("iron",ThreeVector(3*scale,-1*scale,1*scale),ThreeVector(1*scale,1*scale,1*scale)));
  // geometry.AddVolume(SimpleBox("iron",ThreeVector(2*scale,0,0),ThreeVector(1*scale,1*scale,1*scale)));
  geometry.AddVolume(SimpleBox("iron",ThreeVector(5*scale+100,0,0),ThreeVector(10*scale,2*scale,2*scale)));
  
  // Simulator
  Simulator simulator(&geometry);
  Particle electron("electron",11,kElectronMass);
  Bremsstrahlung bremsstrahlung(&electron);
  BetheEnergyLoss bethe_energy_loss;
  PairProduction pairproduction;
  Particle photon("photon",22,0);
  electron.RegisterProcess(&bremsstrahlung);
  // electron.RegisterProcess(&bethe_energy_loss);
  photon.RegisterProcess(&pairproduction);
  simulator.AddParticle(electron);
  simulator.AddParticle(photon);
  simulator.step_size = 0.01;
  simulator.sorting = PARTICLE;
  simulator.thread_multiplier = 1;

  // Arguments
  simulator.device = (runcpu) ? CPU : GPU;
  simulator.debug = true;
  simulator.cpu_threads = 8;
  simulator.pool_size = 200;

  // Tracks
  std::vector<Track> tracks;
  const Float arc = kPi/4;
  const Float E = 10e3; // 100 GeV
  const Float gamma = E / kElectronMass;
  const Float beta = Beta(gamma);
  const Float momentum = sqrt(E*E - kElectronMass*kElectronMass);
  TRandom3 rng((size_t)clock());
  for (int i=0;i<n_tracks;i++) {
    Float angle_phi;
    Float angle_theta;
    if (n_tracks > 1) {
      angle_phi = -arc + 2*arc * (rng.Rndm());
      angle_theta = 2.0 * kPi * (rng.Rndm());
    } else {
      angle_phi = 0;
      angle_theta = 0;
    }
    ThreeVector direction = SphericalToCartesian(momentum,angle_phi,angle_theta);;
    tracks.push_back(Track(11,-1,FourVector(),FourVector(direction[2],direction[0],direction[1],E)));
  }

  // Propagate
  for (int i=0;i<runs;i++) {
    simulator.AddTracks(tracks);
    simulator.steps_per_launch = steps[i];
    simulator.Propagate();
    simulator.ClearTracks();
    benchmark[i] = simulator.GetBenchmark();
  }
  simulator.AddTracks(tracks);
  simulator.device = CPU;
  simulator.Propagate();
  benchmark[runs] = simulator.GetBenchmark();
  std::ofstream outfile;
  outfile.open("Data/bremsstrahlung_benchmark");
  for (int i=0;i<runs;i++) {
    outfile << steps[i] << "," << benchmark[i] << std::endl;
  }
  outfile << 0 << "," << benchmark[runs] << std::endl;
  outfile.close();

  return 0;
}