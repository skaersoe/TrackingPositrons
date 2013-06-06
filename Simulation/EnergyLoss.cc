#include "Geometry/Library.hh"
#include "Simulation/Track.hh"
#include "Simulation/Landau.hh"
#include "Simulation/BetheEnergyLoss.hh"
#include "Simulation/GetTime.hh"

// ROOT crap
#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>
#include "TRandom2.h"
#include "TH1.h"
#include "TH2.h"

using namespace na63;

/**
 * Generates beta*gamma values in a range, generates the parameters for a Landau
 * distribution and samples an amount of values to fill into a histogram.
 */
int main(int argc, char* argv[]) {
  
  // Initialize ROOT crap
  TApplication app("app",&argc,argv);
  TGraph graph;
  TH2F hist_2d(
    "Landau/Vamilov Energy Loss 2D","Bethe/Landau Energy Loss;#beta#gamma;f(dE,dx)",
    500, 0, 100, 500, 0, 500
  );
  // TH1F hist1d("meh","t", 300,0,20);

  const Float kMuonMass = 1.05658372e2; // MeV
  const Float kMuonCharge = -1.0;
  const Float kIronAtomicNumber = 26.0;
  // http://www.physics.nist.gov/cgi-bin/Star/compos.pl?refer=ap&matno=026
  const Float kIronMeanExcitationPotential = 286.0e-6; // eV
  const Float kIronDensity = 7.874; // g/cm^-3
  const Float kIronAtomicWeight = 55.84;
  const Float kCopperAtomicNumber = 29.0;
  const Float kCopperDensity = 8.96;
  const Float kCopperAtomicWeight = 63.546;
  const Float kCopperMeanExcitationPotential = 322.0e-6;

  const unsigned samples_per_energy = 200;

  int i = 0;
  int j = 0;

  TRandom3 rng;
  rng.SetSeed(GetTime().tv_nsec);

  // Parameters to use for run
  const Float mass = kMuonMass;
  const Float charge = kMuonCharge;
  const Float atomic_number = kCopperAtomicNumber;
  const Float density = kCopperDensity;
  const Float atomic_weight = kCopperAtomicWeight;
  const Float mean_excitation_potential = kCopperMeanExcitationPotential;
  const Float dl = 1.0; // cm

  Float gamma = 0.1;
  do {

    i++;
    gamma += 1e-5 * i;

    Float beta = Beta(gamma);
    Float betagamma = beta * gamma;

    // Limits
    if (betagamma != betagamma || betagamma < 0.1) continue;
    if (betagamma > 100.0) break;

    Float dEdx = Bethe_dEdx(beta,mass,charge,atomic_number,
        mean_excitation_potential,density,atomic_weight,dl);
    LandauParameters p = LandauEnergyLossParameters(beta,
        mass,atomic_number,mean_excitation_potential,
        dl);

    //Float kinetic_energy = gamma*kElectronMass - kElectronMass;
    //LandauParameters p = GetElectronCollisionLoss(kIronDensity,kIronAtomicNumber,beta,kinetic_energy,kIronMeanExcitationPotential,-1,dl);
    // LandauParameters p = GetBetheLandauParameters(beta,kMuonMass,kMuonCharge,kIronAtomicNumber,kIronMeanExcitationPotential,dl);

    // Fill graph with mean
    //printf("%f, %f\n",betagamma,dEdx);

    graph.SetPoint(j,betagamma,dl*dEdx);

    // Fill histogram with Landau samples
    for (int i=0;i<samples_per_energy;i++) {
      Float random = ThrowLandauHost(p.mpv,4*p.xi,rng.Rndm());
      hist_2d.Fill(betagamma,random);
    }

    j++;


  } while (1);

  // Save and draw histogram
  // TCanvas canvas_1d;
  // hist_1d.Draw("COLZ");
  //TCanvas canvas_graph;
  TCanvas canvas_2d;
  canvas_2d.SetLogz();
  hist_2d.SetStats(0);
  hist_2d.Draw("COLZ");
  graph.Draw("same l");
  canvas_2d.SaveAs("Plots/stoppingpower_muon.png");

  //TCanvas cccs;
  //hist1d.Draw("hist");
  gPad->WaitPrimitive();

  return 0;
}