#include <iostream>
#include <fstream>
#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <cassert>
#include <sstream>
#include <TLegend.h>

#include "Geometry/Geometry.hh"
#include "Geometry/SimpleBox.hh"
#include "Simulation/GEANT4Bremsstrahlung.hh"
#include "Simulation/GEANT4PairProduction.hh"
#include "Simulation/BetheEnergyLoss.hh"
#include "Simulation/Simulator.hh"
#include "Simulation/GetTime.hh"
#include "Simulation/EnergyLoss.hh"

using namespace na63;

int main(int argc,char *argv[]) {

  const unsigned n_tracks = 64;

  bool plot = false;

  // Run on CPU
  bool runcpu = false;
  for (int i=1;i<argc;i++) {
    std::string arg(argv[i]);
    if (arg == "CPU") {
      runcpu = true;
    }
    if (arg == "plot") {
      plot = true;
    }
  }

  // Geometry
  Geometry geometry;
  const Float scale = 1e3;
  geometry.AddMaterial(Material("vacuum",0.0,0.0,0.0,0.0,0.0));
  geometry.AddMaterial(Material("iron",kIronAtomicNumber,kIronDensity,kIronAtomicWeight,kIronMeanExcitationPotential,kIronRadiationLength));
  geometry.SetBounds(SimpleBox("vacuum",ThreeVector(5*scale,0,0),ThreeVector(10*scale,2*scale,2*scale)));
  // geometry.AddVolume(SimpleBox("iron",ThreeVector(50,1e2+5,0),ThreeVector(30,2e2,4e2)));
  // geometry.AddVolume(SimpleBox("iron",ThreeVector(50,-(1e2+5),0),ThreeVector(30,2e2,4e2)));
  // geometry.AddVolume(SimpleBox("iron",ThreeVector(2*scale,0,0),ThreeVector(1*scale,1*scale,1*scale)));
  // geometry.AddVolume(SimpleBox("iron",ThreeVector(1*scale,1*scale,1*scale),ThreeVector(1*scale,1*scale,1*scale)));
  // geometry.AddVolume(SimpleBox("iron",ThreeVector(3*scale,-1*scale,1*scale),ThreeVector(1*scale,1*scale,1*scale)));
  // geometry.AddVolume(SimpleBox("iron",ThreeVector(2*scale,0,0),ThreeVector(1*scale,1*scale,1*scale)));
  geometry.AddVolume(SimpleBox("iron",ThreeVector(5*scale,0,0),ThreeVector(10*scale,2*scale,2*scale)));

  // Simulator
  Simulator simulator(&geometry);
  Particle electron("electron",11,kElectronMass);
  Bremsstrahlung bremsstrahlung(&electron);
  BetheEnergyLoss bethe_energy_loss;
  PairProduction pairproduction;
  Particle photon("photon",22,0);
  electron.RegisterProcess(&bremsstrahlung);
  electron.RegisterProcess(&bethe_energy_loss);
  photon.RegisterProcess(&pairproduction);
  simulator.AddParticle(electron);
  simulator.AddParticle(photon);
  simulator.step_size = 0.01;
  simulator.sorting = PARTICLE;
  simulator.thread_multiplier = 1;
  const Float xmin = 0;
  const Float xmax = 100;
  const Float xy_dim = 0.4;
  const unsigned xy_bins = 300;
  EnergyLoss energyloss(
    600,xmin,xmax,
    xy_bins,-xy_dim,xy_dim,
    xy_bins,-xy_dim,xy_dim);
  simulator.RecordEnergyLoss(&energyloss);

  // Arguments
  simulator.device = (runcpu) ? CPU : GPU;
  simulator.debug = true;
  simulator.cpu_threads = 8;
  simulator.pool_size = 1000;
  simulator.steps_per_launch = 500;

  // Tracks
  std::vector<Track> tracks;
  const Float arc = kPi/16;
  const Float E = 5e3; // MeV
  const Float gamma = E / kElectronMass;
  simulator.secondary_threshold = 1.0;
  bremsstrahlung.SetSecondaryThreshold(1.0);
  const Float beta = Beta(gamma);
  const Float momentum = sqrt(E*E - kElectronMass*kElectronMass);
  TRandom3 rng((size_t)clock());
  for (int i=0;i<n_tracks;i++) {
    Float angle_phi = 0;
    Float angle_theta = 0.5 * kPi;
    if (n_tracks > 1) {
      // angle_phi = -arc + 2*arc * (rng.Rndm());
      // angle_theta = 2.0 * kPi * (rng.Rndm());
    } else {
      angle_phi = 0.00;
      angle_theta = 0.00;
    }
    ThreeVector direction = SphericalToCartesian(momentum,angle_phi,angle_theta);
    direction.Rotate(ThreeVector(1,0,0));
    tracks.push_back(Track(11,-1,FourVector(),FourVector(direction,E)));
  }
  simulator.AddTracks(tracks);

  // Propagate
  simulator.Propagate();
  unsigned electrons = 0, electrons_live = 0, photons = 0, photons_live = 0;
  Float final_energy = 0;
  tracks = simulator.GetTracks();
  // simulator.PrintTracks();
  // std::ofstream outfile;
  // outfile.open("Data/bremsstrahlung_shower");
  for (int i=0;i<tracks.size();i++) {
    Track t = tracks[i];
    if (t.particle_id == 11) {
      electrons++;
      if (t.momentum[3] > kElectronMass) {
        final_energy += t.momentum[3];
        electrons_live++;
      }
    } else if (t.particle_id == 22) {
      photons++;
      if (t.momentum[3] > 0) {
        final_energy += t.momentum[3];
        photons_live++;
      }
    }
    FourVector vertex = t.vertex();
    // outfile << t.particle_id << ","
    //         << vertex[0] << ","
    //         << vertex[1] << ","
    //         << vertex[2] << ","
    //         << vertex[3] << ","
    //         << t.position[0] << ","
    //         << t.position[1] << ","
    //         << t.position[2] << ","
    //         << t.position[3] << std::endl;
  }
  // outfile.close();
  std::cout << "Simulator returned " << tracks.size() << " particles, of which " << electrons << " are electrons (" << electrons_live << " alive), and " << photons << " are photons (" << photons_live << " alive)." << std::endl;
  Float initial_energy = (Float)n_tracks * E;
  std::cout << "(" << final_energy << " MeV) / (" << initial_energy << " MeV) = " << 100.0*(final_energy/initial_energy) << " percent energy remaining in system." << std::endl;

  if (plot) {
    TApplication app("app",&argc,argv);

    // TH3F *loss_3d = energyloss.Build3DHistogram();
    // loss_3d->SaveAs("Data/energyloss_3d.pdf");

    TH1F* loss_x = energyloss.BuildXHistogram();

    TH2F* loss_xy = energyloss.BuildXYHistogram();
    TH2F* loss_yz = energyloss.BuildYZHistogram();
    loss_x->SetTitle("Energy loss with bremsstrahlung in x;x [cm];Energy [MeV]");
    loss_xy->SetTitle("Energy loss with bremsstrahlung in iron;x [cm];y [cm];Energy [MeV]");
    loss_yz->SetTitle("Energy loss with bremsstrahlung in iron;y [cm];z [cm];Energy [MeV]");
    loss_xy->SetStats(0);
    loss_xy->GetXaxis()->SetLabelSize(0.04);
    loss_xy->GetYaxis()->SetLabelSize(0.04);
    loss_xy->GetZaxis()->SetLabelSize(0.04);
    loss_yz->SetStats(0);
    loss_yz->GetXaxis()->SetLabelSize(0.04);
    loss_yz->GetYaxis()->SetLabelSize(0.04);
    loss_yz->GetZaxis()->SetLabelSize(0.04);
    loss_x->SetStats(0);

    TCanvas *canvas = new TCanvas();
    TLegend *legend = new TLegend(0.63,0.89,0.87,0.80);
    std::stringstream ss;
    ss << "<x> = " << loss_x->GetMean();
    legend->AddEntry("Mean x",ss.str().c_str(),"");
    legend->SetTextSize(0.04);
    canvas->SetLogz();
    canvas->SetRightMargin(0.125);
    loss_xy->Draw("colz");
    legend->Draw();
    // gPad->WaitPrimitive();
    canvas->SaveAs("Plots/energyloss_xy.pdf");

    TCanvas *canvas2 = new TCanvas();
    canvas2->SetLogz();
    canvas2->SetRightMargin(0.125);
    loss_yz->Draw("colz");
    // gPad->WaitPrimitive();
    canvas2->SaveAs("Plots/energyloss_yz.pdf");
  }

  // Analyze
  // TApplication app("app",&argc,argv);
  // TCanvas canvas;
  // TH1F electrons("Bremsstrahlung electrons","Bremsstrahlung electrons",96,-arc,arc);
  // TH1F photons("Bremsstrahlung photons","Bremsstrahlung photons",96,-arc,arc);
  // electrons.SetLineColor(kRed);
  // photons.SetLineColor(kBlue);
  // std::ofstream outfile;
  // outfile.open("Data/bremsstrahlung_shower");
  // FourVector zero_vertex;
  // int i=0;
  // while (i<simulator.TrackSize()) {
  //   Track t = simulator.GetTrack(i);
  //   FourVector vertex = t.vertex();
  //   // if (t.particle_id == 11 && vertex == zero_vertex) {
  //   //   std::cout << i << ": " << t.particle_id << ", " << t.position << ", " << t.momentum << std::endl;
  //   // }
  //   if (t.particle_id == 11) {
  //     electrons.Fill(angles[t.initial_index]);
  //   } else if (t.particle_id == 22) {
  //     photons.Fill(angles[t.initial_index]);
  //   }
  //   // outfile << t.particle_id << ","
  //   //         << vertex[0] << ","
  //   //         << vertex[1] << ","
  //   //         << vertex[2] << ","
  //   //         << vertex[3] << ","
  //   //         << t.position[0] << ","
  //   //         << t.position[1] << ","
  //   //         << t.position[2] << ","
  //   //         << t.position[3] << std::endl;
  //   i++;
  // }
  // outfile.close();
  // std::cout << "Copied back " << i << " tracks." << std::endl;
  // photons.Draw();
  // electrons.Draw("same");
  // photons.SetTitle("Bremsstrahlung");
  // photons.GetXaxis()->SetTitle("Angle [radians]");
  // canvas.Modified();
  // canvas.Update();
  //canvas.SaveAs("Plots/bremsstrahlung.png");
  //gPad->WaitPrimitive();

  return 0;
}