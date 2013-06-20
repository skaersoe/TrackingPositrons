#include <iostream>
#include <cmath>
#include <TApplication.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TH1F.h>
#include <TF1.h>
#include <TStyle.h>
#include <TRoot.h>
#include "Simulation/GetTime.hh"

int main(int argc, char* argv[]) {

  TH1F *range = new TH1F(
    "Interaction probability for dx=0.1, X_0=10",
    "Interaction probability;Distance travelled before interaction [cm]",
    200,
    0,
    100
  );

  TRandom3 rng(na63::InSeconds(na63::GetTime()));

  float dx = 0.05;
  const float X_0 = 10;
  const unsigned runs = 1e6;

  const float probability = 1 - exp(-dx/X_0);

  for (int i=0;i<runs;i++) {
    int j = 0;
    do {
      j++;
    } while (rng.Rndm() > probability);
    range->Fill(j * dx);
  }
  TApplication app("app",&argc,argv);
  gROOT->SetBatch(kTRUE);
  TCanvas *canvas = new TCanvas("","",1280,1024);
  range->Draw();
  range->Fit("expo");
  TF1 *f = range->GetFunction("expo");
  f->SetLineColor(kBlue);
  // f->SetLineStyle(7);
  f->Draw("same");
  gStyle->SetOptFit();
  gPad->WaitPrimitive();
  printf("Saving...\n");
  canvas->SaveAs("Plots/probability.pdf");
  return 0;
}