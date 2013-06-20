#include <iostream>
#include <fstream>
#include <TApplication.h>
#include <TCanvas.h>
#include "Simulation/Simulator.hh"
#include "Geometry/Geometry.hh"
#include "Geometry/SimpleBox.hh"
#include "Simulation/BetheEnergyLoss.hh"
#include "Simulation/GEANT4Bremsstrahlung.hh"
#include "Simulation/GEANT4PairProduction.hh"
#include "Simulation/EnergyLoss.hh"

using namespace na63;

double RunSimulator(Simulator *simulator) {
  simulator->Propagate();
  simulator->ClearTracks();
  return simulator->GetBenchmark();
}

int main(int argc,char** argv) {

  // Benchmark vector
  std::vector<double> benchmark;

  // Set up geometry
  Geometry geometry;
  Float scale = 1e3;
  geometry.AddMaterial(Material("vacuum",0.0,0.0,0.0,0.0,0.0));
  geometry.AddMaterial(Material("iron",kIronAtomicNumber,kIronDensity,kIronAtomicWeight,kIronMeanExcitationPotential,kIronRadiationLength));
  geometry.SetBounds(SimpleBox("vacuum",ThreeVector(5*scale,0,0),ThreeVector(10*scale,2*scale,2*scale)));
  geometry.AddVolume(SimpleBox("iron",ThreeVector(5*scale+100,0,0),ThreeVector(10*scale,2*scale,2*scale)));

  // Set up simulator
  Simulator simulator(&geometry);
  simulator.debug = false;
  simulator.step_size = 0.01;
  simulator.pool_size = 0;
  simulator.sorting = PARTICLE;
  simulator.cpu_threads = 8;
  simulator.steps_per_launch = 500;
  simulator.record_energyloss = false;
  Particle electron("electron",11,kElectronMass);
  Particle photon("photon",22,0);
  Bremsstrahlung bremsstrahlung(&electron);
  BetheEnergyLoss bethe_energy_loss;
  PairProduction pairproduction;
  electron.RegisterProcess(&bethe_energy_loss);

  // Get command line input
  std::string which_benchmark("");
  for (int i=1;i<argc;i++) {
    which_benchmark += argv[i];
  }
  std::cout << "Received input \"" << which_benchmark << "\"..." << std::endl;

  // Run appropriate benchmark
  if (which_benchmark.find("energy") != std::string::npos && which_benchmark.find("bremsstrahlung") != std::string::npos) {
    std::ofstream outfile;
    electron.RegisterProcess(&bremsstrahlung);
    photon.RegisterProcess(&pairproduction);
    simulator.AddParticle(electron);
    simulator.AddParticle(photon);
    simulator.gpu_bremsstrahlung = true;
    simulator.thread_multiplier = 128;
    outfile.open("Data/benchmark_energy_bremsstrahlung");
    const int runs = 9;
    const int n_tracks = 32;
    const Float energies[11] = {1e1,5e1,1e2,5e2,1e3,5e3,1e4,5e4,1e5,5e5,1e6};
    outfile << "type,energy,runtime,heap size" << std::endl;
    for (int r=0;r<runs;r++) {
      std::vector<Track> tracks;
      const Float energy = energies[r];
      simulator.pool_size = 0.15 * energy;
      const Float momentum = sqrt(energy*energy - kElectronMass*kElectronMass);
      for (int i=0;i<n_tracks;i++) {
        tracks.push_back(Track(11,-1,FourVector(),FourVector(momentum,0,0,energy)));
      }
      std::cout << "Running test " << r+1 << "/" << runs << " for energy " << energy << "..." << std::endl;
      simulator.device = CPU;
      for (int i=0;i<3;i++) {
        if (i > 1) simulator.pool_size = 0.15 * energies[runs-1];
        simulator.AddTracks(tracks);
        double runtime = RunSimulator(&simulator);
        benchmark.push_back(runtime);
        outfile << i << "," << energy << "," << runtime << "," << simulator.pool_size << std::endl;
        simulator.device = GPU;
      }
    }
    outfile.close();
    return 0;
  }

  if (which_benchmark.find("energy") != std::string::npos) {
    simulator.AddParticle(electron);
    const int n_tracks = 4096;
    std::ofstream outfile;
    outfile.open("Data/benchmark_energy");
    const int runs = 9;
    const Float energies[] = {1e1,5e1,1e2,5e2,1e3,5e3,1e4,5e4,1e5,5e5,1e6};
    simulator.gpu_bremsstrahlung = false;
    outfile << "energy,runtime" << std::endl;
    for (int r=0;r<runs;r++) {
      std::vector<Track> tracks;
      const Float energy = energies[r];
      const Float momentum = sqrt(energy*energy - kElectronMass*kElectronMass);
      for (int i=0;i<n_tracks;i++) {
        tracks.push_back(Track(11,-1,FourVector(),FourVector(momentum,0,0,energy)));
      }
      std::cout << "Running test " << r+1 << "/" << runs << " for energy " << energy << "..." << std::endl;
      simulator.device = CPU;
      simulator.step_size = 0.01;
      simulator.steps_per_launch = 5000;
      for (int i=0;i<2;i++) {
        simulator.AddTracks(tracks);
        double runtime = RunSimulator(&simulator);
        benchmark.push_back(runtime);
        outfile << i << "," << energy << "," << runtime << std::endl;
        simulator.device = GPU;
      }
    }
    outfile.close();
    return 0;
  }

  if (which_benchmark.find("steps") != std::string::npos) {
    int n_tracks = 4096;
    std::ofstream outfile;
    if (which_benchmark.find("bremsstrahlung") != std::string::npos) {
      outfile.open("Data/benchmark_steps_bremsstrahlung");
      n_tracks = 64;
      electron.RegisterProcess(&bremsstrahlung);
      photon.RegisterProcess(&pairproduction);
      simulator.AddParticle(photon);
      simulator.gpu_bremsstrahlung = true;
      simulator.thread_multiplier = 8;
    } else {
      outfile.open("Data/benchmark_steps");
      simulator.gpu_bremsstrahlung = false;
    }
    simulator.AddParticle(electron);
    const int runs = 8;
    const Float steps[] = {100,200,500,1000,2000,3000,4000,5000};
    outfile << "steps per launch,runtime" << std::endl;
    const Float energy = 1e3;
    simulator.pool_size = 0.15 * energy;
    const Float momentum = sqrt(energy*energy - kElectronMass*kElectronMass);
    std::vector<Track> tracks;
    simulator.device = GPU;
    simulator.step_size = 0.1;
    for (int i=0;i<n_tracks;i++) {
      tracks.push_back(Track(11,-1,FourVector(),FourVector(momentum,0,0,energy)));
    }
    double runtime;
    for (int r=0;r<runs;r++) {
      std::cout << "Running test " << r+1 << "/" << runs << " for " << steps[r] << " steps per launch..." << std::endl;
      simulator.steps_per_launch = steps[r];
      simulator.AddTracks(tracks);
      runtime = RunSimulator(&simulator);
      outfile << steps[r] << "," << runtime << std::endl;
    }
    simulator.AddTracks(tracks);
    simulator.device = CPU;
    runtime = RunSimulator(&simulator);
    outfile << "0" << "," << runtime << std::endl;
    outfile.close();
    return 0;
  }

  if (which_benchmark.find("incident") != std::string::npos) {
    const Float energy = 10e3;
    const Float momentum = energy*energy - kElectronMass*kElectronMass;
    const Track track(11,-1,FourVector(),FourVector(momentum,0,0,energy));
    std::vector<Track> tracks;
    const int runs = 9;
    const int n_tracks[] = {64,128,256,512,1024,2048,4096,8192,16384};
    // electron.RegisterProcess(&bremsstrahlung);
    // photon.RegisterProcess(&pairproduction);
    // simulator.AddParticle(photon);
    simulator.AddParticle(electron);
    simulator.gpu_bremsstrahlung = false;
    std::ofstream outfile;
    outfile.open("Data/benchmark_incident");
    for (int r=0;r<runs;r++) {
      std::cout << "Running test " << r+1 << "/" << runs << " for " << n_tracks[r] << " incident particles." << std::endl;
      // simulator.pool_size = (int)(0.1250 * energy * (Float)n_tracks[r]);
      // simulator.thread_multiplier = 4096 / n_tracks[r];
      // if (simulator.thread_multiplier > 16) simulator.thread_multiplier = 16;
      for (int i=0;i<n_tracks[r];i++) {
        tracks.push_back(track);
      }
      simulator.device = CPU;
      for (int i=0;i<2;i++) {
        double runtime;
        simulator.AddTracks(tracks);
        runtime = RunSimulator(&simulator);
        outfile << i << "," << n_tracks[r] << "," << runtime << std::endl;
        simulator.device = GPU;
      }
      tracks.clear();
    }
    outfile.close();
    return 0;
  }

  if (which_benchmark.find("step_size") != std::string::npos) {
    // electron.RegisterProcess(&bremsstrahlung);
    // photon.RegisterProcess(&pairproduction);
    // simulator.AddParticle(photon);
    Particle muon("muon",13,kMuonMass);
    muon.RegisterProcess(&bethe_energy_loss);
    simulator.AddParticle(muon);
    Geometry geometry2;
    scale = 1e3;
    geometry2.AddMaterial(Material("vacuum",0.0,0.0,0.0,0.0,0.0));
    geometry2.AddMaterial(Material("copper",kCopperAtomicNumber,kCopperDensity,kCopperAtomicWeight,kCopperMeanExcitationPotential,kCopperRadiationLength));
    geometry2.AddMaterial(Material("iron",kIronAtomicNumber,kIronDensity,kIronAtomicWeight,kIronMeanExcitationPotential,kIronRadiationLength));
    geometry2.AddMaterial(Material("lead",kLeadAtomicNumber,kLeadDensity,kLeadAtomicWeight,kLeadMeanExcitationPotential,kLeadRadiationLength));
    geometry2.AddMaterial(Material("graphite",kGraphiteAtomicNumber,kGraphiteDensity,kGraphiteAtomicWeight,kGraphiteMeanExcitationPotential,kGraphiteRadiationLength));
    simulator.geometry = &geometry2;
    // EnergyLoss energyloss(
    //   800,0,80,
    //   300,-1.0,1.0,
    //   300,-1.0,1.0
    // );
    // simulator.RecordEnergyLoss(&energyloss);
    simulator.debug = false;
    simulator.pool_size = 0;
    simulator.thread_multiplier = 1;
    simulator.steps_per_launch = 2500;
    const int runs = 6;
    const Float step_size[] = {0.2,0.1,0.05,0.025,0.01,0.005};
    const Float energy = 10e3;
    const Float momentum = energy*energy - kMuonMass*kMuonMass;
    const int n_tracks = 2048;
    Track track(13,-1,FourVector(),FourVector(momentum,0,0,energy));
    std::vector<Track> tracks;
    for (int i=0;i<n_tracks;i++) {
      tracks.push_back(track);
    }
    std::ofstream outfile;
    outfile.open("Data/benchmark_stepsize");
    for (int i=0;i<10;i++) {
      geometry2.SetBounds(SimpleBox("vacuum",ThreeVector(),ThreeVector(2.0*(scale*(i+1)),10*scale,10*scale)));
      geometry2.AddVolume(SimpleBox("lead",ThreeVector(scale*(i + 0.125),0,0),ThreeVector(0.25*scale,10*scale,10*scale)));
      geometry2.AddVolume(SimpleBox("graphite",ThreeVector(scale*(i + 0.375),0,0),ThreeVector(0.25*scale,10*scale,10*scale)));
      geometry2.AddVolume(SimpleBox("copper",ThreeVector(scale*(i + 0.625),0,0),ThreeVector(0.25*scale,10*scale,10*scale)));
      geometry2.AddVolume(SimpleBox("iron",ThreeVector(scale*(i + 0.875),0,0),ThreeVector(0.25*scale,10*scale,10*scale)));
    }
    for (int r=0;r<runs;r++) {
      std::cout << "Running with step size " << step_size[r] << "..." << std::endl;
      simulator.step_size = step_size[r];
      simulator.device = CPU;
      for (int i=0;i<2;i++) {
        simulator.AddTracks(tracks);
        double runtime = RunSimulator(&simulator);
        outfile << i << "," << step_size[r] << "," << runtime << std::endl;
        simulator.device = GPU;
      }
    }
    outfile.close();
    // TH2F *hist = energyloss.BuildXYHistogram();
    // TApplication app("app",&argc,argv);
    // TCanvas canvas;
    // canvas.SetLogz();
    // hist->Draw("colz");
    // gPad->WaitPrimitive();
    return 0;
  }

  if (which_benchmark.find("threshold") != std::string::npos) {
    electron.RegisterProcess(&bremsstrahlung);
    photon.RegisterProcess(&pairproduction);
    simulator.AddParticle(photon);
    simulator.AddParticle(electron);
    const Float energy = 5e3;
    const Float momentum = energy*energy - kElectronMass*kElectronMass;
    const Track track(11,-1,FourVector(),FourVector(momentum,0,0,energy));
    std::vector<Track> tracks;
    const int n_tracks = 128;
    simulator.thread_multiplier = 32;
    for (int i=0;i<n_tracks;i++) {
      tracks.push_back(track);
    }
    const int runs = 7;
    const Float threshold[] = {1.0,2.0,5.0,10.0,50.0,100.0,energy};
    const int heap_size[] = {800,700,650,550,500,500,500};
    std::ofstream outfile;
    outfile.open("Data/benchmark_threshold");
    simulator.debug = true;
    for (int r=0;r<runs;r++) {
      std::cout << "Running for threshold " << threshold[r] << "..." << std::endl;
      bremsstrahlung.SetSecondaryThreshold(2.0 * kElectronMass * threshold[r]);
      simulator.secondary_threshold = 2.0 * kElectronMass * threshold[r];
      simulator.device = GPU;
      simulator.pool_size = heap_size[r];
      for (int i=0;i<2;i++) {
        simulator.AddTracks(tracks);
        double runtime = RunSimulator(&simulator);
        outfile << i << "," << threshold[r] << "," << runtime << std::endl;
        simulator.device = CPU;
      }
    }
    outfile.close();
    return 0;
  }

  std::cerr << "Invalid input. No benchmark performed." << std::endl;
  return -1; 
}