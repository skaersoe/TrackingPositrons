#include <iostream>
#include <fstream>
#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <TRandom3.h>
#include <cassert>

#include "Geometry/Geometry.hh"
#include "Geometry/Box.hh"
#include "Simulation/Bremsstrahlung.hh"

using namespace na63;

int main(int argc,char *argv[]) {

  const unsigned n_tracks = 256;

  // Geometry
  Geometry geometry;
  geometry.AddMaterial(Material("vacuum",0.0,0.0,0.0));
  geometry.AddMaterial(Material("iron",26.0,286.0,13.84));
  geometry.SetBounds(Box("vacuum",ThreeVector(2e2,0,0),ThreeVector(4e2,4e2,4e2)));
  geometry.AddVolume(Box("iron",ThreeVector(2e2,0,0),ThreeVector(1e2,1e2,1e2)));

  // Simulator
  Simulator simulator(&geometry);
  Particle electron("electron",11,kElectronMass);
  Particle photon("photon",22,0);
  electron.RegisterProcess(Bremsstrahlung);
  photon.RegisterProcess(Bremsstrahlung);
  simulator.AddParticle(electron);
  simulator.AddParticle(photon);
  simulator.step_size = 0.1;
  simulator.sorting = X;

  // Arguments
  simulator.device = CPU;
  simulator.debug = true;
  simulator.cpu_threads = 8;

  // Tracks
  std::vector<Track> tracks;
  const Float arc = kPi/4;
  const Float E = 10; // 100 MeV
  const Float gamma = E / kElectronMass;
  const Float beta = Beta(gamma);
  Float angles[n_tracks];
  TRandom3 rng((size_t)clock());
  for (int i=0;i<n_tracks;i++) {
    Float angle;
    if (n_tracks > 1) {
      angle = -arc + 2*arc * (rng.Rndm());
    } else {
      angle = 0;
    }
    angles[i] = angle;
    Float bx = beta * cos(angle);
    Float by = beta * sin(angle);
    Float px = Gamma(bx) * kElectronMass * bx;
    Float py = Gamma(by) * kElectronMass * by;
    Float pz = 0;
    tracks.push_back(Track(11,-1,FourVector(),FourVector(px,py,pz,E)));
  }
  simulator.AddTracks(tracks);

  // Propagate
  simulator.Propagate();

  // Analyze
  TApplication app("app",&argc,argv);
  TCanvas canvas;
  TH1F electrons("Bremsstrahlung electrons","Bremsstrahlung electrons",96,-arc,arc);
  TH1F photons("Bremsstrahlung photons","Bremsstrahlung photons",96,-arc,arc);
  electrons.SetLineColor(kRed);
  photons.SetLineColor(kBlue);
  std::ofstream outfile;
  outfile.open("Data/bremsstrahlung_shower");
  for (int i=0;i<simulator.TrackSize();i++) {
    Track t = simulator.GetTrack(i);
    if (t.particle_id == 11) {
      electrons.Fill(angles[t.initial_index]);
    } else if (t.particle_id == 22) {
      photons.Fill(angles[t.initial_index]);
    }
    FourVector vertex = t.vertex();
    outfile << t.particle_id << ","
            << vertex[0] << ","
            << vertex[1] << ","
            << vertex[2] << ","
            << vertex[3] << ","
            << t.position[0] << ","
            << t.position[1] << ","
            << t.position[2] << ","
            << t.position[3] << std::endl;
  }
  outfile.close();
  photons.Draw();
  electrons.Draw("same");
  photons.SetTitle("Bremsstrahlung");
  photons.GetXaxis()->SetTitle("Angle [radians]");
  canvas.Modified();
  canvas.Update();
  //canvas.SaveAs("Plots/bremsstrahlung.png");
  //gPad->WaitPrimitive();

  return 0;
}