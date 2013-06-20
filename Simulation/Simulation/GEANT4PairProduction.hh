#ifndef NA63_SIMULATION_GEANT4PAIRPRODUCTION_H
#define NA63_SIMULATION_GEANT4PAIRPRODUCTION_H

#include "Simulation/Process.hh"

namespace na63 {
  
class PairProduction : public Process {

public:
  PairProduction();
  void SampleSecondaries(Track* track, const Material* couple);
  virtual void Query(Track* track, const Material* material, const Float dl);

private:
  void CalcLPMFunctions(Float k, Float eplusEnergy);
  void SetupForMaterial(const Material* material_new);

  const Material *material;
  Float gLPM;
  Float phiLPM;
  Float xiLPM;
  Float lpmEnergy;
  bool fLPMflag;

  static const Float xgi[8], wgi[8];
  static const Float Fel_light[5];
  static const Float Finel_light[5];
  static const Float facFel;
  static const Float facFinel;

  static const Float preS1, logTwo;

};

} // End namespace na63

#endif /* NA63_SIMULATION_GEANT4PAIRPRODUCTION_H */