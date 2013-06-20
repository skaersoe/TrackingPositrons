#include <cassert>
#include <iostream>
#include <vector>
#include <cmath>
#include <thread>
#include <sstream>
#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>
#include "TRandom2.h"
#include "TH2.h"

using namespace std;

typedef float Float;

typedef struct {
  Float mpv;
  Float xi;
} LandauParameters;

LandauParameters LandauEnergyLossParameters(const Float beta,
    const Float mass, const Float atomic_number,
    const Float mean_excitation_potential, const Float dl) {

  LandauParameters p;

  // Calculate necessary values
  Float beta_squared = beta*beta;
  Float gamma = 1/sqrt(1-beta_squared);
  Float gamma_squared = gamma*gamma;
  p.xi = 0.5 * 0.307075 * atomic_number * dl / beta_squared;

  p.mpv = p.xi * (log(2.0 * mass * beta_squared
      * gamma_squared / mean_excitation_potential)
      + log(p.xi/mean_excitation_potential) + 0.200 - beta_squared);

  return p;

}

int main(int argc, char* argv[]) {

  TApplication app("app",&argc,argv);
  TGraph *graph_electron = new TGraph();
  TGraph *graph_muon = new TGraph();

  Float kElectronMass = 0.5;
  Float kMuonMass = 100;

  for (int i=1;i<1000;i++) {
    Float energy = (Float)i;
    Float gamma_e = (kElectronMass + energy) / kElectronMass;
    Float gamma_m = (kMuonMass + energy) / kMuonMass;
    Float beta_e = sqrt(1-1/pow(gamma_e,2));
    Float beta_m = sqrt(1-1/pow(gamma_m,2));
    assert(beta_e == beta_e);
    assert(beta_m == beta_m);
    LandauParameters p_e = LandauEnergyLossParameters(beta_e,kElectronMass,26.0,0.000286,0.1);
    LandauParameters p_m = LandauEnergyLossParameters(beta_m,kMuonMass,26.0,0.000286,0.1);
    assert(p_e.mpv == p_e.mpv);
    assert(p_m.mpv == p_m.mpv);
    assert(energy > 0);
    graph_electron->SetPoint(i,energy,p_e.mpv/energy);
    graph_muon->SetPoint(i,energy,p_m.mpv/energy);
    if (i < 900) continue;
    printf("Electron lost %f / %f energy\n",p_e.mpv,energy);
    printf("Muon lost %f / %f energy\n",p_m.mpv,energy);
  }

  graph_electron->Draw("al");
  graph_muon->Draw("same l");
  graph_electron->SetLineColor(kBlue);
  gPad->WaitPrimitive();

  return 0;
}