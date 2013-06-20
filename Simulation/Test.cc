#include <iostream>
#include "Geometry/Constants.hh"
#include "Simulation/GEANT4Bremsstrahlung.hh"

using namespace na63;
using namespace std;

void KillTrack(Track& track) {
  track.Kill();
}

int main(void) {
  Geometry *geometry = new Geometry();
  Material *copper = new Material("copper",kCopperAtomicNumber,kCopperDensity,kCopperMeanExcitationPotential,kCopperRadiationLength);
  geometry->AddMaterial(*copper);
  Material *lead = new Material("lead",kLeadAtomicNumber,kLeadDensity,kLeadMeanExcitationPotential,kLeadRadiationLength);
  geometry->AddMaterial(*lead);
  Material *iron = new Material("iron",kIronAtomicNumber,kIronDensity,kIronMeanExcitationPotential,kIronRadiationLength);
  geometry->AddMaterial(*iron);
  Simulator *simulator = new Simulator(geometry);
  Particle *electron = new Particle("electron",11,kElectronMass);
  simulator->AddParticle(*electron);
  Particle *photon = new Particle("photon",22,0);
  simulator->AddParticle(*photon);
  Float energy = 100;
  Float momentum = sqrt(energy*energy - electron->mass()*electron->mass());
  simulator->AddTrack(Track(11,-1,FourVector(),FourVector(momentum,0,0,energy)));
  Bremsstrahlung b(electron);
  while (simulator->TrackSize() == 1) {
    b.Query(simulator->GetTrackAddress(0),copper,0.1);
  }
  simulator->PrintTracks();
  return 0;
}