#include <iostream>
#include <vector>

// ROOT crap
#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>
#include "TRandom2.h"
#include "TH1.h"
#include "TH2.h"

#include "Simulation/Simulator.hh"
#include "Geometry/Geometry.hh"
#include "Geometry/Box.hh"
#include "Geometry/Sphere.hh"
#include "Simulation/GetTime.hh"

#include "Simulation/BetheEnergyLoss.hh"

using namespace na63;

int main(int argc,char *argv[]) {

  // Set up geometry
  Geometry geometry;
  geometry.AddMaterial(Material("vacuum",0,0));
  geometry.AddMaterial(Material("iron",26,286));
  geometry.SetBounds(Box("vacuum",ThreeVector(0,0,0),ThreeVector(1e10,1e10,1e10)));
  geometry.AddVolume(Box("iron",ThreeVector(0,0,0),ThreeVector(1e10,1e10,1e10)));

  // Create simulator
  Simulator simulator = Simulator(&geometry);
  Particle muon = Particle("muon",13,-1,105.6583715);
  muon.RegisterProcess(BetheEnergyLoss);
  simulator.AddParticle(muon);
  simulator.step_size = 0.1;

  // Set parameters
  simulator.device = CPU;
  simulator.debug = false;
  simulator.cpu_threads = 8;
  const unsigned x_samples = 256;
  const unsigned y_samples = 4;
  const unsigned n_tracks = x_samples * y_samples;
  const Float energy_range[2] = {1e3,1e5};

  // Generate tracks
  std::vector<Track> tracks;
  Float energy[x_samples];
  for (int i=0;i<x_samples;i++) {
    energy[i] = energy_range[0] + (energy_range[1] - energy_range[0]) * ((Float)i / (Float)x_samples);
    for (int j=0;j<y_samples;j++) {
      tracks.push_back(Track(13,FourVector(),
          FourVector(energy[i]/kC,0,0,energy[i]/kC)));
    }
  }
  simulator.AddTracks(tracks);

  // Propagate
  simulator.Propagate();

  // Report
  TApplication app("app",&argc,argv);
  TCanvas canvas;
  TGraph graph;
  /*TH2F hist_2d(
      "Penetration length","Penetration length;E [GeV];x [cm]",
      100,energy_range[0],energy_range[1],100,0,500
  );*/
  for (int i=0;i<x_samples;i++) {
    Track t = simulator.GetTrack(i);
    graph.SetPoint(i,1e-3*energy[i],t.position[0]);
    //hist_2d.Fill(energy[i],t.position[0]);
  }
  //hist_2d.Draw("COLZ");
  graph.Draw("same ap");
  graph.SetMarkerStyle(5);
  graph.SetTitle("Muon depth in iron");
  graph.GetXaxis()->SetTitle("Energy [GeV]");
  graph.GetYaxis()->SetTitle("Depth [cm]");
  graph.GetXaxis()->CenterTitle();
  graph.GetYaxis()->CenterTitle();
  gPad->WaitPrimitive();

  return 0;

}