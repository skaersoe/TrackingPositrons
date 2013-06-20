#include <iostream>
#include <vector>
#include <sstream>

// ROOT crap
#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>
#include "TRandom3.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"

#include "Simulation/Simulator.hh"
#include "Geometry/Geometry.hh"
#include "Geometry/Box.hh"
#include "Geometry/Sphere.hh"
#include "Simulation/GetTime.hh"

#include "Simulation/BetheEnergyLoss.hh"

using namespace na63;

int main(int argc,char *argv[]) {

  TRandom3 rng(InSeconds(GetTime()));

  // Set up geometry
  Geometry geometry;
  geometry.AddMaterial(Material("vacuum",0,0,0,0,0));
  geometry.AddMaterial(Material("iron",kIronAtomicNumber,kIronDensity,kIronAtomicWeight,kIronMeanExcitationPotential,kIronRadiationLength));
  geometry.SetBounds(Sphere("vacuum",ThreeVector(0,0,0),1e10));
  geometry.AddVolume(Sphere("iron",ThreeVector(0,0,0),1e10));

  // Create simulator
  Simulator simulator = Simulator(&geometry);
  Particle muon = Particle("muon",13,kMuonMass);
  BetheEnergyLoss bethe;
  muon.RegisterProcess(&bethe);
  simulator.AddParticle(muon);
  simulator.step_size = 0.1;

  // Set parameters
  simulator.device = CPU;
  simulator.debug = true;
  simulator.cpu_threads = 8;
  simulator.steps_per_launch = 1000;
  simulator.step_size = 0.01;
  const unsigned x_samples = 128;
  const unsigned y_samples = 8;
  const unsigned n_tracks = x_samples * y_samples;
  const Float energy_range[2] = {1e3,1e5};
  const Float energy_step = (energy_range[1] - energy_range[0]) / x_samples;

  // Generate tracks
  std::vector<Track> tracks;
  Float energy[n_tracks];
  for (int i=0;i<x_samples;i++) {
    for (int j=0;j<y_samples;j++) {
      int idx = i*y_samples+j;
      energy[idx] = (Float)i * energy_step + rng.Rndm() * energy_step;
      tracks.push_back(Track(13,-1,FourVector(),
          FourVector(energy[idx]/kMuonMass,0,0,energy[idx])));
    }
  }

  // Propagate on CPU
  simulator.AddTracks(tracks);
  double elapsed_cpu = InSeconds(GetTime());
  simulator.Propagate();
  elapsed_cpu = InSeconds(GetTime()) - elapsed_cpu;
  std::vector<Track> result_cpu = simulator.GetTracks();
  simulator.ClearTracks();
  // Propagate on GPU
  simulator.device = GPU;
  simulator.AddTracks(tracks);
  double elapsed_gpu = InSeconds(GetTime());
  simulator.Propagate();
  elapsed_gpu = InSeconds(GetTime()) - elapsed_gpu;
  std::vector<Track> result_gpu = simulator.GetTracks();

  // Report
  TApplication app("app",&argc,argv);
  TCanvas canvas;
  TGraph graph_cpu;
  TGraph graph_gpu;
  /*TH2F hist_2d(
      "Penetration length","Penetration length;E [GeV];x [cm]",
      100,energy_range[0],energy_range[1],100,0,500
  );*/
  for (int i=0;i<n_tracks;i++) {
    graph_gpu.SetPoint(i,1e-3*energy[result_gpu[i].initial_index],result_gpu[i].position[0]);
    graph_cpu.SetPoint(i,1e-3*energy[i],result_cpu[i].position[0]);
    //hist_2d.Fill(energy[i],t.position[0]);
  }
  //hist_2d.Draw("COLZ");
  graph_cpu.Draw("AP");
  graph_gpu.Draw("same P");
  graph_cpu.SetMarkerStyle(1);
  graph_gpu.SetMarkerStyle(1);
  graph_cpu.SetMarkerColor(4);
  graph_gpu.SetMarkerColor(8);
  graph_cpu.SetTitle("Muon depth in iron");
  graph_cpu.GetXaxis()->SetTitle("Energy [GeV]");
  graph_cpu.GetYaxis()->SetTitle("Depth [cm]");
  graph_cpu.GetXaxis()->CenterTitle();
  graph_cpu.GetYaxis()->CenterTitle();
  canvas.Modified();
  canvas.Update();
  TLegend *legend = new TLegend(0.13,0.89,0.33,0.79);
  std::stringstream cpu_legend;
  cpu_legend << "CPU, ran in " << elapsed_cpu << " seconds.";
  std::stringstream gpu_legend;
  gpu_legend << "GPU, ran in " << elapsed_gpu << " seconds.";
  legend->AddEntry(&graph_cpu,cpu_legend.str().c_str(),"p");
  legend->AddEntry(&graph_gpu,gpu_legend.str().c_str(),"p");
  legend->Draw();
  gPad->WaitPrimitive();
  //canvas.SaveAs("Plots/depth_cpu_gpu.png");

  return 0;

}