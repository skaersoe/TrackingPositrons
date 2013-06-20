#include <iostream>
#include <fstream>
#include "Simulation/Simulator.hh"
#include "Geometry/Geometry.hh"
#include "Geometry/SimpleBox.hh"
#include "Simulation/GEANT4Bremsstrahlung.hh"
#include "Simulation/BetheEnergyLoss.hh"
#include "Simulation/GEANT4PairProduction.hh"

using namespace std;
using namespace na63;

int main(int argc, char** argv) {

  const int n_tracks = 32;
  const int runs = 11;
  const Float energy[] = {1e1,5e1,1e2,5e2,1e3,5e3,1e4,5e4,1e5,5e5,1e6};
  const Float mass = kElectronMass;
  int *result = new int[runs];

  Geometry geometry;
  geometry.AddMaterial(Material("vacuum",0,0,0,0,0));
  geometry.AddMaterial(Material("iron",kIronAtomicNumber,kIronDensity,kIronAtomicWeight,kIronMeanExcitationPotential,kIronRadiationLength));
  geometry.SetBounds(SimpleBox("vacuum",ThreeVector(0,0,0),ThreeVector(1e4,1e4,1e4)));
  geometry.AddVolume(SimpleBox("iron",ThreeVector(0,0,0),ThreeVector(1e4,1e4,1e4)));

  Simulator simulator(&geometry);
  Particle electron("electron",11,kElectronMass);
  Bremsstrahlung bremsstrahlung(&electron);
  BetheEnergyLoss bethe_energy_loss;
  PairProduction pairproduction;
  Particle photon("photon",22,0);
  electron.RegisterProcess(&bremsstrahlung);
  electron.RegisterProcess(&bethe_energy_loss);
  photon.RegisterProcess(&pairproduction);
  simulator.AddParticle(electron);
  simulator.AddParticle(photon);
  simulator.device = CPU;
  simulator.debug = false;
  simulator.cpu_threads = 8;

  // std::ofstream outfile;
  // outfile.open("Data/particles_per_energy");
  // for (int r=0;r<runs;r++) {

  //   std::cout << "Propagating with energy " << energy[r] << " MeV..." << std::endl;

  //   Float momentum = sqrt(energy[r]*energy[r] - mass*mass);
  //   Track track(11,-1,FourVector(),FourVector(momentum,0,0,energy[r]));

  //   for (int p=0;p<n_tracks;p++) {

  //     simulator.AddTrack(track);

  //     simulator.Propagate();

  //     outfile << energy[r] << "," << simulator.TrackSize() << std::endl;

  //     simulator.ClearTracks();

  //   }
  // }
  // outfile.close();

  std::ofstream threshold_outfile;
  threshold_outfile.open("Data/particles_per_threshold");
  const int threshold_runs = 8;
  const Float threshold_energy = 5e3;
  const Float threshold[] = {1.0,2.0,5.0,10.0,25.0,50.0,100.0,threshold_energy};
  const Float threshold_momentum = sqrt(threshold_energy*threshold_energy - kElectronMass*kElectronMass);
  const Track threshold_track(11,-1,FourVector(),FourVector(threshold_momentum,0,0,threshold_energy));
  for (int r=0;r<threshold_runs;r++) {
    std::cout << "Propagating with threshold " << threshold[r] << "..." << std::endl;
    bremsstrahlung.SetSecondaryThreshold(threshold[r]);
    for (int i=0;i<n_tracks;i++) {
      simulator.AddTrack(threshold_track);
      simulator.Propagate();
      threshold_outfile << threshold[r] << "," << simulator.TrackSize() << std::endl;
      simulator.ClearTracks();
    }
  }
  threshold_outfile.close();

  return 0;
}