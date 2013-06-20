#include <cmath>
#include <algorithm>
#include <TRandom3.h>

#include "Geometry/Constants.hh"
#include "Simulation/GetTime.hh"
#include "Simulation/GEANT4PairProduction.hh"
#include "Simulation/ScreenFunctions.hh"
#include "Simulation/Track.hh"
#include "Geometry/Material.hh"

namespace na63 {

PairProduction::PairProduction() {
  fLPMflag = true;
  material = nullptr;
  gLPM = 0;
  phiLPM = 0;
  xiLPM = 0;
  lpmEnergy = 0;
}

TRandom3 rng_pairproduction(InSeconds(GetTime()));

const Float PairProduction::facFel = log(184.15);
const Float PairProduction::facFinel = log(1194.); // 1440.

const Float PairProduction::preS1 = 1./(184.15*184.15);
const Float PairProduction::logTwo = log(2.);

const Float PairProduction::xgi[]={ 0.0199, 0.1017, 0.2372, 0.4083,
             0.5917, 0.7628, 0.8983, 0.9801 };
const Float PairProduction::wgi[]={ 0.0506, 0.1112, 0.1569, 0.1813,
             0.1813, 0.1569, 0.1112, 0.0506 };
const Float PairProduction::Fel_light[]  = {0., 5.31  , 4.79  , 4.74 ,  4.71};
const Float PairProduction::Finel_light[] = {0., 6.144 , 5.621 , 5.805 , 5.924};

inline Float Phi1(Float delta) {
   Float screenVal;

   if (delta > 1.)
     screenVal = 21.12 - 4.184*std::log(delta+0.952);
   else
     screenVal = 20.868 - delta*(3.242 - 0.625*delta);

   return screenVal;
}

inline Float Phi2(Float delta) {
   Float screenVal;

   if (delta > 1.)
     screenVal = 21.12 - 4.184*std::log(delta+0.952);
   else
     screenVal = 20.209 - delta*(1.930 + 0.086*delta);

   return screenVal;
}

void PairProduction::SetupForMaterial(const Material* material_new) {
  if (material == material_new) return;
  material = material_new;
  lpmEnergy = material->radiation_length()*material->density()*kLPMConstant;
}

void PairProduction::SampleSecondaries(
    Track* track,
    const Material* couple) {

// The secondaries e+e- energies are sampled using the Bethe - Heitler
// cross sections with Coulomb correction.
// A modified version of the random number techniques of Butcher & Messel
// is used (Nuc Phys 20(1960),15).
//
// GEANT4 internal units.
//
// Note 1 : Effects due to the breakdown of the Born approximation at
//          low energy are ignored.
// Note 2 : The differential cross section implicitly takes account of 
//          pair creation in both nuclear and atomic electron fields.
//          However triplet prodution is not generated.

  Float GammaEnergy = track->kinetic_energy();
  ThreeVector GammaDirection = track->momentum.Normalized();

  Float epsil;
  Float epsil0 = kElectronMass/GammaEnergy;
  if (epsil0 > 1.0) return;

  SetupForMaterial(couple);

  assert(material != nullptr);

  // do it fast if GammaEnergy < 2. MeV
  static const Float Egsmall = 2.0 * MeV;

  // select randomly one element constituing the material
  // const G4Element* anElement = SelectRandomAtom(aMaterial, theGamma, GammaEnergy);

  if (GammaEnergy < Egsmall) {

    epsil = epsil0 + (0.5-epsil0)*rng_pairproduction.Rndm();

  } else {
    // now comes the case with GammaEnergy >= 2. MeV

    // Extract Coulomb factor for this Element
    Float Z3 = std::pow(couple->atomic_number(),3);
    Float FZ = 8.*(Z3);
    if (GammaEnergy > 50.*MeV) FZ += 8.*(couple->coulomb_correction());

    // limits of the screening variable
    Float screenfac = 136.*epsil0/(Z3);
    Float screenmax = std::exp((42.24 - FZ)/8.368) - 0.952 ;
    Float temp = 4.0*screenfac;
    
    Float screenmin = std::min(temp,screenmax);

    // limits of the energy sampling
    Float epsil1 = 0.5 - 0.5*sqrt(1. - screenmin/screenmax) ;
    Float epsilmin = std::max(epsil0,epsil1) , epsilrange = 0.5 - epsilmin;

    //
    // sample the energy rate of the created electron (or positron)
    //
    //Float epsil, screenvar, greject ;
    Float  screenvar, greject;

    Float F10 = ScreenFunction1(screenmin) - FZ;
    Float F20 = ScreenFunction2(screenmin) - FZ;
    temp = 0.0;
    Float NormF1 = std::max(F10*epsilrange*epsilrange,temp); 
    Float NormF2 = std::max(1.5*F20,0.);

    do {
      if (NormF1/(NormF1+NormF2) > rng_pairproduction.Rndm()) {
        epsil = 0.5 - epsilrange*pow(rng_pairproduction.Rndm(), 0.333333);
        screenvar = screenfac/(epsil*(1-epsil));
        if (fLPMflag && GammaEnergy>100.*GeV) {
          CalcLPMFunctions(GammaEnergy,GammaEnergy*epsil);
          greject = xiLPM*((gLPM+2.*phiLPM)*Phi1(screenvar) - gLPM*Phi2(screenvar) - phiLPM*FZ)/F10;
        } else {
          greject = (ScreenFunction1(screenvar) - FZ)/F10;
        }
                    
      } else { 
        epsil = epsilmin + epsilrange*rng_pairproduction.Rndm();
        screenvar = screenfac/(epsil*(1-epsil));
        if (fLPMflag && GammaEnergy>100.*GeV) {
          CalcLPMFunctions(GammaEnergy,GammaEnergy*epsil);
          greject = xiLPM*((0.5*gLPM+phiLPM)*Phi1(screenvar) + 0.5*gLPM*Phi2(screenvar) - 0.5*(gLPM+phiLPM)*FZ)/F20;
          // printf("gLPM: %f, phiLPM: %f, xiLPM: %f\n",gLPM,phiLPM,xiLPM);
        } else {
          greject = (ScreenFunction2(screenvar) - FZ)/F20;
        }
      }

    } while (greject < rng_pairproduction.Rndm());

  }   //  end of epsil sampling
   
  //
  // fixe charges randomly
  //

  Float ElectTotEnergy, PositTotEnergy;
  if (rng_pairproduction.Rndm() > 0.5) {

    ElectTotEnergy = (1.-epsil)*GammaEnergy;
    PositTotEnergy = epsil*GammaEnergy;
     
  } else {
    
    PositTotEnergy = (1.-epsil)*GammaEnergy;
    ElectTotEnergy = epsil*GammaEnergy;
  }

  //
  // scattered electron (positron) angles. ( Z - axis along the parent photon)
  //
  //  universal distribution suggested by L. Urban 
  // (Geant3 manual (1993) Phys211),
  //  derived from Tsai distribution (Rev Mod Phys 49,421(1977))

  Float u;
  const Float a1 = 0.625 , a2 = 3.*a1 , d = 27. ;

  if (9./(9.+d) > rng_pairproduction.Rndm()) u= - log(rng_pairproduction.Rndm()*rng_pairproduction.Rndm())/a1;
  else u= - log(rng_pairproduction.Rndm()*rng_pairproduction.Rndm())/a2;

  Float TetEl = u*kElectronMass/ElectTotEnergy;
  Float TetPo = u*kElectronMass/PositTotEnergy;
  Float Phi  = 2.0 * kPi * rng_pairproduction.Rndm();
  Float dxEl= sin(TetEl)*cos(Phi),dyEl= sin(TetEl)*sin(Phi),dzEl=cos(TetEl);
  Float dxPo=-sin(TetPo)*cos(Phi),dyPo=-sin(TetPo)*sin(Phi),dzPo=cos(TetPo);
   
  //
  // kinematic of the created pair
  //
  // the electron and positron are assumed to have a symetric
  // angular distribution with respect to the Z axis along the parent photon.

  Float temp = 0.0;
  Float ElectKineEnergy = std::max(temp,ElectTotEnergy - kElectronMass);

  ThreeVector ElectDirection(dxEl, dyEl, dzEl);
  ElectDirection.Rotate(GammaDirection);
  ElectDirection.Extend(sqrt(ElectTotEnergy*ElectTotEnergy - kElectronMass*kElectronMass));

  // create G4DynamicParticle object for the particle1  
  // G4DynamicParticle* aParticle1= new G4DynamicParticle(
  //        theElectron,ElectDirection,ElectKineEnergy);
  Track aParticle1(11,-1,track->position,FourVector(ElectDirection,ElectTotEnergy));
  
  // the e+ is always created (even with Ekine=0) for further annihilation.

  temp = 0.0;
  Float PositKineEnergy = std::max(temp,PositTotEnergy - kElectronMass);

  ThreeVector PositDirection(dxPo, dyPo, dzPo);
  PositDirection.Rotate(GammaDirection);  
  PositDirection.Extend(sqrt(PositTotEnergy*PositTotEnergy - kElectronMass*kElectronMass));;


  // create G4DynamicParticle object for the particle2 
  // G4DynamicParticle* aParticle2= new G4DynamicParticle(
  //                     thePositron,PositDirection,PositKineEnergy);
  Track aParticle2(11,1,track->position,FourVector(PositDirection,PositTotEnergy));

  // Kill photon
  track->Kill();

  // Spawn children
  std::vector<Track> children;
  children.push_back(aParticle1);
  children.push_back(aParticle2);
  assert((abs(children[0].momentum[3] + children[1].momentum[3] - GammaEnergy)) < 1e-4);
  track->SpawnChildren(children);
}

void PairProduction::CalcLPMFunctions(Float k, Float eplusEnergy) {


  // *** calculate lpm variable s & sprime ***
  // Klein eqs. (78) & (79)
  Float sprime = sqrt(0.125*k*lpmEnergy/(eplusEnergy*(k-eplusEnergy)));

  Float Z = material->atomic_number();
  Float s1 = preS1*pow(Z,0.6666666);
  Float logS1 = 2./3.*log(Z)-2.*facFel;
  Float logTS1 = kLogTwo+logS1;

  xiLPM = 2.;

  if (sprime>1) 
    xiLPM = 1.;
  else if (sprime>sqrt(2.)*s1) {
    Float h  = log(sprime)/logTS1;
    xiLPM = 1+h-0.08*(1-h)*(1-sqrt(1-h))/logTS1;
  }

  Float s0 = sprime/sqrt(xiLPM); 
  //   G4cout<<"k="<<k<<" y="<<eplusEnergy/k<<G4endl;
  //   G4cout<<"s0="<<s0<<G4endl;
  
  // *** calculate supression functions phi and G ***
  // Klein eqs. (77)
  Float s2=s0*s0;
  Float s3=s0*s2;
  Float s4=s2*s2;

  if (s0<0.1) {
    // high suppression limit
    phiLPM = 6.*s0 - 18.84955592153876*s2 + 39.47841760435743*s3 
      - 57.69873135166053*s4;
    gLPM = 37.69911184307752*s2 - 236.8705056261446*s3 + 807.7822389*s4;
  }
  else if (s0<1.9516) {
    // intermediate suppression
    // using eq.77 approxim. valid s0<2.      
    phiLPM = 1.-exp(-6.*s0*(1.+(3.-kPi)*s0)
    +s3/(0.623+0.795*s0+0.658*s2));
    if (s0<0.415827397755) {
      // using eq.77 approxim. valid 0.07<s<2
      Float psiLPM = 1-exp(-4*s0-8*s2/(1+3.936*s0+4.97*s2-0.05*s3+7.50*s4));
      gLPM = 3*psiLPM-2*phiLPM;
    }
    else {
      // using alternative parametrisiation
      Float pre = -0.16072300849123999 + s0*3.7550300067531581 + s2*-1.7981383069010097 
  + s3*0.67282686077812381 + s4*-0.1207722909879257;
      gLPM = tanh(pre);
    }
  }
  else {
    // low suppression limit valid s>2.
    phiLPM = 1. - 0.0119048/s4;
    gLPM = 1. - 0.0230655/s4;
  }

  // *** make sure suppression is smaller than 1 ***
  // *** caused by Migdal approximation in xi    ***
  if (xiLPM*phiLPM>1. || s0>0.57)  xiLPM=1./phiLPM;

}

void PairProduction::Query(Track* track, const Material* material,
    const Float dl) {

  assert(material != nullptr);

  if (track->energy() < 2.0*kElectronMass) return;


  // Use radiation length for probability
  Float chance_to_interact = 1 - std::exp(-dl/material->radiation_length());
  if (rng_pairproduction.Rndm() > chance_to_interact) return;

  SampleSecondaries(track,material);

}


} // End namespace na63