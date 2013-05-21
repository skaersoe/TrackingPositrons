#include <iostream>
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

  const unsigned n_tracks = 512;

  // Geometry
  Geometry geometry;
  geometry.AddMaterial(Material("vacuum",0.0,0.0));
  geometry.AddMaterial(Material("iron",26.0,286.0));
  geometry.SetBounds(Box("vacuum",ThreeVector(2e3,0,0),ThreeVector(4e3,4e3,4e3)));
  geometry.AddVolume(Box("iron",ThreeVector(2e3,0,0),ThreeVector(1e3,1e3,1e3)));

  // Simulator
  Simulator simulator(&geometry);
  Particle electron = Particle("electron",11,-1,kElectronMass);
  InitializeBremsstrahlung();
  electron.RegisterProcess(Bremsstrahlung);
  simulator.AddParticle(electron);
  simulator.AddParticle(Particle("photon",22,0,0));
  simulator.step_size = 0.1;

  // Arguments
  simulator.device = CPU;
  simulator.debug = true;
  simulator.cpu_threads = 8;

  // Tracks
  std::vector<Track> tracks;
  const Float arc = kPi/4;
  const Float beta = 0.95;
  const Float E = Gamma(beta) * kElectronMass;
  Float angles[n_tracks];
  TRandom3 rng((size_t)clock());
  for (int i=0;i<n_tracks;i++) {
    Float angle = -arc + 2*arc * (rng.Rndm());
    angles[i] = angle;
    Float bx = beta * cos(angle);
    Float by = beta * sin(angle);
    Float px = Gamma(bx) * kElectronMass * bx;
    Float py = Gamma(by) * kElectronMass * by;
    Float pz = 0;
    tracks.push_back(Track(11,FourVector(),FourVector(px,py,pz,E)));
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
  for (int i=0;i<simulator.TrackSize();i++) {
    Track t = simulator.GetTrack(i);
    if (t.particle_id == 11) {
      electrons.Fill(angles[t.initial_index]);
    } else if (t.particle_id == 22) {
      photons.Fill(angles[t.initial_index]);
    }
  }
  photons.Draw();
  electrons.Draw("same");
  photons.SetTitle("Bremsstrahlung");
  photons.GetXaxis()->SetTitle("Angle [radians]");
  canvas.Modified();
  canvas.Update();
  gPad->WaitPrimitive();

  return 0;
}