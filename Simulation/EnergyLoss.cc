#include <cassert>

#include "Geometry/Library.hh"
#include "Simulation/Track.hh"
#include "Simulation/Landau.hh"
#include "Simulation/BetheEnergyLoss.hh"
#include "Simulation/GetTime.hh"
#include "Geometry/Constants.hh"

// ROOT crap
#include <TRoot.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>
#include "TRandom2.h"
#include "TH1.h"
#include "TH2.h"
#include <TAxis.h>
#include <TMath.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TColor.h>

using namespace na63;

TRandom3 rng;

const unsigned samples_per_energy = 1 << 15;

void LogAxisBins(TAxis *axis) {

  int n_bins = axis->GetNbins();

  Float min = axis->GetXmin();
  Float max = axis->GetXmax();

  assert(max > min);

  assert(n_bins > 1);

  Float *bins_new = new Float[n_bins];

  for (int i=0;i<n_bins;i++) {
    bins_new[i] = min + (max - min) * pow((Float)i / (Float)(n_bins - 1),10);
    assert(bins_new[i] >= min);
    if (i > 0) {
      assert(bins_new[i] >= bins_new[i-1]);
    }
  }

  axis->Set(n_bins-1,bins_new);
  //axis->Set(n_bins,axis->GetXbins());

}

void EnergyLoss(TH2F *hist, TGraph *graph, TGraph *graph_mpv, const Float mass, const Float charge,
    const Float atomic_number, const Float mean_excitation_potential,
    const Float density, const Float atomic_weight, const Float dl) {

  TAxis *xaxis = hist->GetXaxis();
  int bins = xaxis->GetNbins();
  for (int i=0;i<bins;i++) {

    Float betagamma = hist->GetBinCenter(i);
    Float gamma = sqrt(betagamma*betagamma + 1);

    Float dEdx = Bethe_dEdx(gamma,mass,charge,atomic_number,
        mean_excitation_potential,density,atomic_weight,dl);
    LandauParameters p = LandauEnergyLossParameters(gamma,
        mass,atomic_number,density,atomic_weight,mean_excitation_potential,
        dl);

    //Float kinetic_energy = gamma*kElectronMass - kElectronMass;
    //LandauParameters p = GetElectronCollisionLoss(kIronDensity,kIronAtomicNumber,beta,kinetic_energy,kIronMeanExcitationPotential,-1,dl);
    // LandauParameters p = GetBetheLandauParameters(beta,kMuonMass,kMuonCharge,kIronAtomicNumber,kIronMeanExcitationPotential,dl);

    // Fill graph with mean
    //printf("%f, %f\n",betagamma,dEdx);

    Float energy = gamma * mass;
    Float momentum = sqrt(energy*energy - mass*mass);

    if (graph != nullptr)
      graph->SetPoint(i,betagamma,dEdx/density);
    if (graph_mpv != nullptr)
      graph_mpv->SetPoint(i,betagamma,p.mpv/(dl*density));

    // Fill histogram with Landau samples
    for (int j=0;j<samples_per_energy;j++) {
      Float random = ThrowLandauHost(p.mpv,4*p.xi,rng.Rndm()) / (dl * density);
      hist->Fill(betagamma,random);
    }

  }
}

/**
 * Generates beta*gamma values in a range, generates the parameters for a Landau
 * distribution and samples an amount of values to fill into a histogram.
 */
int main(int argc, char* argv[]) {

  rng.SetSeed(GetTime().tv_nsec);

  // Parameters to use for run
  const Float mass = kMuonMass;;
  const Float charge = -1;
  const Float atomic_number = kLeadAtomicNumber;
  const Float density = kLeadDensity;
  const Float atomic_weight = kLeadAtomicWeight;
  const Float mean_excitation_potential = kLeadMeanExcitationPotential;
  const Float dl = 0.001; // cm
  
  // MDJ contour
  // const int NRGBs = 7, NCont = 999;
  // gStyle->SetNumberContours(NCont);
  // Double_t stops[NRGBs] = { 0.05, 0.10, 0.25, 0.45, 0.60, 0.75, 1.00 };
  // Double_t red[NRGBs]   = { 1.00, 0.00, 0.00, 0.00, 0.97, 0.97, 0.10 };
  // Double_t green[NRGBs] = { 1.00, 0.97, 0.30, 0.40, 0.97, 0.00, 0.00 };
  // Double_t blue[NRGBs]  = { 1.00, 0.97, 0.97, 0.00, 0.00, 0.00, 0.00 };
  // TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);

  // Initialize ROOT crap
  TApplication app("app",&argc,argv);
  TGraph *graph_muon = new TGraph();
  TGraph *graph_mpv = new TGraph();
  // TGraph *graph_electron= new TGraph();
  TH2F *hist_2d = new TH2F(
    "Muon energy loss in copper","Muon energy loss in copper;#beta#gamma;Stopping power [MeV cm^{2} / g]",
    450, 0.1, 100, 450, 1e-1, 1e2
  );
  // hist_2d->SetMaximum(1.5e3);
  hist_2d->SetContour(256);
  LogAxisBins(hist_2d->GetXaxis());
  LogAxisBins(hist_2d->GetYaxis());
  hist_2d->GetXaxis()->SetRangeUser(0.1*(1.001),1e2);
  // hist_2d->GetYaxis()->SetRangeUser(2e-1,0.9e2);

  EnergyLoss(hist_2d,
             graph_muon,
             graph_mpv,
             mass,
             charge,
             atomic_number,
             mean_excitation_potential,
             density,
             atomic_weight,
             dl);

  // EnergyLoss(hist_2d,
  //            nullptr,
  //            nullptr,
  //            kProtonMass,
  //            kProtonCharge,
  //            atomic_number,
  //            mean_excitation_potential,
  //            density,
  //            atomic_weight,
  //            dl);

  // Save and draw histogram
  // TCanvas canvas_1d;
  // hist_1d.Draw("COLZ");
  //TCanvas canvas_graph;


  TCanvas *canvas = new TCanvas();
  canvas->SetLogx();
  canvas->SetLogy();
  canvas->SetLogz();
  hist_2d->SetStats(0);
  hist_2d->Draw("COLZ");
  graph_muon->Draw("same l");
  graph_muon->SetLineStyle(1);
  graph_muon->SetLineWidth(3);
  graph_mpv->Draw("same l");
  graph_mpv->SetLineStyle(7);
  graph_mpv->SetLineWidth(3);
  graph_mpv->SetLineColor(kBlue);
  // graph_electron->Draw("same l");
  // graph_electron->SetLineStyle(7);
  // graph_electron->SetLineWidth(3);
  TLegend *legend = new TLegend(0.47,0.89,0.89,0.74);
  legend->AddEntry(graph_muon,"Mean energy loss","l");
  legend->AddEntry(graph_mpv,"Most probable energy loss","l");
  legend->AddEntry("Step size","dx = 1e-3cm","");
  legend->SetTextSize(0.04);
  legend->SetFillColor(kWhite);
  legend->Draw();
  canvas->Update();
  gPad->WaitPrimitive();
  printf("Saving...\n");
  canvas->SaveAs("Plots/stoppingpower_muon.pdf");

  //TCanvas cccs;
  //hist1d.Draw("hist");

  return 0;
}