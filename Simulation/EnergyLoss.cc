#include "Geometry/Library.hh"
#include "Simulation/Track.hh"
#include "Simulation/BetheEnergyLoss.hh"

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
    "Bethe/Landau Energy Loss 2D","Bethe/Landau Energy Loss;#beta#gamma;dE/dx",
    500, 0, 100, 500, 0, 100
  );
  TH1F hist1d("meh","t", 300,0,20);

  const Float kMuonMass = 1.05658372e2; // MeV
  const Float kMuonCharge = -1.0;
  const Float kIronAtomicNumber = 26.0;
  // http://www.physics.nist.gov/cgi-bin/Star/compos.pl?refer=ap&matno=026
  const Float kIronMeanExcitationPotential = 286.0e-6; // eV
  const Float dl = 0.001; // cm

  const unsigned samples_per_energy = 200;

  int i = 0;
  int j = 0;
  do {
    i++;

    Float beta = 1e-6 * i;
    Float betagamma = beta * Gamma(beta);

    // Limits
    if (betagamma < 0.1) continue;
    if (betagamma > 100.0) break;
    j++;

    LandauParameters p = GetSkewedLandauParameters(beta,
        kMuonMass,kMuonCharge,kIronAtomicNumber,kIronMeanExcitationPotential,
        dl);

    // Fill graph with mean
    graph.SetPoint(i,betagamma,p.mean);

    // Fill histogram with Landau samples
    for (int i=0;i<samples_per_energy;i++) {
      Float random = rng.Landau(p.mean,p.sigma);
      // hist_1d.Fill(random);
      hist_2d.Fill(betagamma,random);

      if (betagamma > 1 and betagamma < 1.02)
      {
        hist1d.Fill(random);
      }
    }


  } while (1);

  // Save and draw histogram
  // TCanvas canvas_1d;
  // hist_1d.Draw("COLZ");
  //TCanvas canvas_graph;
  TCanvas canvas_2d;
  hist_2d.Draw("COLZ");
  graph.Draw("same l");
  //canvas.SaveAs("Plots/dE_dx.png");

  TCanvas cccs;
  hist1d.Draw("hist");
  gPad->WaitPrimitive();

  return 0;
}