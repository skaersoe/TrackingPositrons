#include <sstream>

#include "Geometry/Geometry.hh"
#include "Geometry/SimpleBox.hh"
#include "Simulation/Simulator.hh"
#include "Simulation/EnergyLoss.hh"
#include "Simulation/GEANT4Bremsstrahlung.hh"
#include "Simulation/GEANT4PairProduction.hh"
#include "Simulation/BetheEnergyLoss.hh"
#include "Simulation/GetTime.hh"

#include <TApplication.h>
#include <TCanvas.h>

using namespace na63;
using namespace std;

int main(int argc, char* argv[]) {
  
  // Set up geometry
  Geometry geometry;
  geometry.AddMaterial(Material("vacuum",0,0,0,0,0));
  Material iron(
    "iron",
    kIronAtomicNumber,
    kIronDensity,
    kIronAtomicWeight,
    kIronMeanExcitationPotential,
    kIronRadiationLength
  );
  Material plastic(
    "plastic", // Polyvinyltoluene
    0.54141,
    1.032,
    // Use weighted mean of constituents
    (1.008 * 0.085 + 12.01 * 0.915) / 2,
    64.7,
    42.54
  );
  geometry.AddMaterial(iron);
  geometry.AddMaterial(plastic);
  SimpleBox bounds(
    "vacuum",
    ThreeVector(1*m,0,0),
    ThreeVector(2*m,1*m,1*m)
  );
  geometry.SetBounds(bounds);
  const Float layer_size = 5 * cm;
  const int layers = 20;
  for (int l=0;l<layers;l++) {
    geometry.AddVolume(
      SimpleBox(
        "iron",
        ThreeVector((2.0*l + 0.5)*layer_size,0,0),
        ThreeVector(layer_size,1*m,1*m)
      )
    );
    geometry.AddVolume(
      SimpleBox(
        "plastic",
        ThreeVector((2.0*l + 1.0 + 0.5)*layer_size,0,0),
        ThreeVector(layer_size,1*m,1*m)
      )
    );
  }

  // Set up simulator
  Simulator simulator(&geometry);
  simulator.device = CPU;
  simulator.cpu_threads = 8;
  simulator.record_energyloss = true;

  // Set up energy loss recording
  EnergyLoss energy_loss(
    600,0,200, // x
    300,-5,5,  // y
    300,-5,5   // z
  );
  simulator.RecordEnergyLoss(&energy_loss);

  // Set up particle types
  Particle electron("electron",11,kElectronMass);
  Particle photon("photon",22,0.0);
  Bremsstrahlung bremsstrahlung(&electron);
  PairProduction pair_production;
  BetheEnergyLoss bethe_energy_loss;
  electron.RegisterProcess(&bremsstrahlung);
  electron.RegisterProcess(&bethe_energy_loss);
  photon.RegisterProcess(&pair_production);
  simulator.AddParticle(electron);
  simulator.AddParticle(photon);

  // Create initial tracks
  const int n_tracks = 64;
  const Float energy = 50 * GeV;
  const Float mass = kElectronMass;
  const Float momentum = sqrt(energy*energy - mass*mass);
  vector<Track> tracks(
    // Create 64 identical electrons
    64,
    Track(
      11,
      -1,
      FourVector(),
      FourVector(momentum,0,0,energy)
    )
  );
  simulator.AddTracks(tracks);

  // Run propagation, get benchmark
  simulator.Propagate();
  double benchmark = simulator.GetBenchmark();

  // Plot xy-plane of energy loss and save to file
  TCanvas c;
  c.SetRightMargin(0.125);
  c.SetLogz();
  TH2F *hist_xy = energy_loss.BuildXYHistogram();
  stringstream title;
  title << "Bremsstrahlung in iron/plastic layers, "
        << "propagated in " << benchmark
        << " seconds.;x [cm];y [cm];Energy [MeV]";
  hist_xy->SetTitle(title.str().c_str());
  hist_xy->SetStats(0);
  hist_xy->GetXaxis()->SetRangeUser(0,200);
  hist_xy->Draw("colz");
  c.SaveAs("Plots/layer_energyloss.pdf");

  return 0;
}