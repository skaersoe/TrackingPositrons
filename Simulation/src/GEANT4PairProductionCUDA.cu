#include "Geometry/LibraryCUDA.cuh"
#include "Geometry/Constants.hh"
#include "Simulation/ScreenFunctionsCUDA.cuh"
#include "Simulation/TrackGPU.cuh"
#include "Simulation/Spawner.cuh"
#include <cfloat>

namespace na63 {

namespace {

__device__ __constant__ Float facFel = 5.21575064289;

__device__ __constant__ Float preS1 = 1./(184.15*184.15);

__device__ __constant__ Float Egsmall = 2.0 * MeV;

__device__ inline
Float Phi1(const Float delta) {
   Float screenVal;

   if (delta > 1.)
     screenVal = 21.12 - 4.184*logf(delta+0.952);
   else
     screenVal = 20.868 - delta*(3.242 - 0.625*delta);

   return screenVal;
}

__device__ inline
Float Phi2(const Float delta) {
   Float screenVal;

   if (delta > 1.)
     screenVal = 21.12 - 4.184*logf(delta+0.952);
   else
     screenVal = 20.209 - delta*(1.930 + 0.086*delta);

   return screenVal;
}

typedef struct {
  Float xiLPM;
  Float phiLPM;
  Float gLPM;
} LPMFunctions;

__device__
LPMFunctions CalcLPMFunctions(
    const Float k,
    const Float eplusEnergy,
    const Float lpmEnergy,
    const Float Z) {

  LPMFunctions ret;

  // *** calculate lpm variable s & sprime ***
  // Klein eqs. (78) & (79)
  Float sprime = sqrt(0.125*k*lpmEnergy/(eplusEnergy*(k-eplusEnergy)));

  Float s1 = preS1*pow(Z,Float(0.6666666));
  Float logS1 = 2./3.*logf(Z)-2.*facFel;
  Float logTS1 = 0.69314718056 + logS1;

  ret.xiLPM = 2.;

  if (sprime>1) 
    ret.xiLPM = 1.;
  else if (sprime>sqrt(2.)*s1) {
    Float h  = logf(sprime)/logTS1;
    ret.xiLPM = 1+h-0.08*(1-h)*(1-sqrt(1-h))/logTS1;
  }

  Float s0 = sprime/sqrt(ret.xiLPM); 
  //   G4cout<<"k="<<k<<" y="<<eplusEnergy/k<<G4endl;
  //   G4cout<<"s0="<<s0<<G4endl;
  
  // *** calculate supression functions phi and G ***
  // Klein eqs. (77)
  Float s2=s0*s0;
  Float s3=s0*s2;
  Float s4=s2*s2;

  if (s0<0.1) {
    // high suppression limit
    ret.phiLPM = 6.*s0 - 18.84955592153876*s2 + 39.47841760435743*s3 
      - 57.69873135166053*s4;
    ret.gLPM = 37.69911184307752*s2 - 236.8705056261446*s3 + 807.7822389*s4;
  }
  else if (s0<1.9516) {
    // intermediate suppression
    // using eq.77 approxim. valid s0<2.      
    ret.phiLPM = 1.-exp(-6.*s0*(1.+(3.-kPi)*s0)
    +s3/(0.623+0.795*s0+0.658*s2));
    if (s0<0.415827397755) {
      // using eq.77 approxim. valid 0.07<s<2
      Float psiLPM = 1-exp(-4*s0-8*s2/(1+3.936*s0+4.97*s2-0.05*s3+7.50*s4));
      ret.gLPM = 3*psiLPM-2*ret.phiLPM;
    }
    else {
      // using alternative parametrisiation
      Float pre = -0.16072300849123999 + s0*3.7550300067531581 + s2*-1.7981383069010097 
  + s3*0.67282686077812381 + s4*-0.1207722909879257;
      ret.gLPM = tanh(pre);
    }
  }
  else {
    // low suppression limit valid s>2.
    ret.phiLPM = 1. - 0.0119048/s4;
    ret.gLPM = 1. - 0.0230655/s4;
  }

  // *** make sure suppression is smaller than 1 ***
  // *** caused by Migdal approximation in xi    ***
  if (ret.xiLPM*ret.phiLPM>1. || s0>0.57) ret.xiLPM=1./ret.phiLPM;

  return ret;

}

__device__
void SampleSecondaries(
    GPUTrack* track,
    // MDJ: replaces MaterialCutsCouple kinda wrong, be aware...
    const ParticlePars* particle,
    const MaterialPars* material,
    const int index,
    const int child_index,
    curandState *rng_state,
    Float cut_energy = 2*kElectronMass,
    Float max_energy = FLT_MAX,
    bool fLPMflag = true) {

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

  Float GammaEnergy = track->momentum[3] - particle->mass;
  GPUThreeVector GammaDirection;
  ThreeVector_Normalized(GammaDirection,track->momentum);

  Float epsil;
  Float epsil0 = kElectronMass/GammaEnergy;
  if (epsil0 > 1.0) return;

  Float lpmEnergy = material->radiation_length*material->density*kLPMConstant;

  // select randomly one element constituing the material
  // const G4Element* anElement = SelectRandomAtom(aMaterial, theGamma, GammaEnergy);

  if (GammaEnergy < Egsmall) {

    epsil = epsil0 + (0.5-epsil0)*curand_uniform(rng_state);

  } else {
    // now comes the case with GammaEnergy >= 2. MeV

    // Extract Coulomb factor for this Element
    Float Z3 = pow(material->atomic_number,3);
    Float FZ = 8.*(Z3);
    if (GammaEnergy > 50.*MeV) FZ += 8.*(material->coulomb_correction);

    // limits of the screening variable
    Float screenfac = 136.*epsil0/(Z3);
    Float screenmax = exp((42.24 - FZ)/8.368) - 0.952 ;
    Float temp = 4.0*screenfac;
    
    Float screenmin = min(temp,screenmax);

    // limits of the energy sampling
    Float epsil1 = 0.5 - 0.5*sqrt(1. - screenmin/screenmax) ;
    Float epsilmin = max(epsil0,epsil1) , epsilrange = 0.5 - epsilmin;

    //
    // sample the energy rate of the created electron (or positron)
    //
    //Float epsil, screenvar, greject ;
    Float  screenvar, greject;

    Float F10 = CUDA_ScreenFunction1(screenmin) - FZ;
    Float F20 = CUDA_ScreenFunction2(screenmin) - FZ;
    temp = 0.0;
    Float NormF1 = max(F10*epsilrange*epsilrange,temp); 
    Float NormF2 = max(1.5*F20,0.);

    LPMFunctions lpm;

    do {
      if (NormF1/(NormF1+NormF2) > curand_uniform(rng_state)) {
        epsil = 0.5 - epsilrange*pow(curand_uniform(rng_state),Float(0.333333));
        screenvar = screenfac/(epsil*(1-epsil));
        if (fLPMflag && GammaEnergy>100.*GeV) {
          lpm = CalcLPMFunctions(GammaEnergy,GammaEnergy*epsil,lpmEnergy,Z);
          greject = lpm.xiLPM*((lpm.gLPM+2.*lpm.phiLPM)*Phi1(screenvar) - lpm.gLPM*Phi2(screenvar) - lpm.phiLPM*FZ)/F10;
        } else {
          greject = (CUDA_ScreenFunction1(screenvar) - FZ)/F10;
        }
                    
      } else { 
        epsil = epsilmin + epsilrange*curand_uniform(rng_state);
        screenvar = screenfac/(epsil*(1-epsil));
        if (fLPMflag && GammaEnergy>100.*GeV) {
          lpm = CalcLPMFunctions(GammaEnergy,GammaEnergy*epsil,lpmEnergy,Z);
          greject = lpm.xiLPM*((0.5*lpm.gLPM+lpm.phiLPM)*Phi1(screenvar) + 0.5*lpm.gLPM*Phi2(screenvar) - 0.5*(lpm.gLPM+lpm.phiLPM)*FZ)/F20;
          // printf("gLPM: %f, phiLPM: %f, xiLPM: %f\n",gLPM,phiLPM,xiLPM);
        } else {
          greject = (CUDA_ScreenFunction2(screenvar) - FZ)/F20;
        }
      }

    } while (greject < curand_uniform(rng_state));

  }   //  end of epsil sampling
   
  //
  // fixe charges randomly
  //

  Float ElectTotEnergy, PositTotEnergy;
  if (curand_uniform(rng_state) > 0.5) {

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

  if (9./(9.+d) > curand_uniform(rng_state)) u= - logf(curand_uniform(rng_state)*curand_uniform(rng_state))/a1;
  else u= - logf(curand_uniform(rng_state)*curand_uniform(rng_state))/a2;

  Float TetEl = u*kElectronMass/ElectTotEnergy;
  Float TetPo = u*kElectronMass/PositTotEnergy;
  Float Phi  = 2.0 * kPi * curand_uniform(rng_state);
  Float dxEl= sin(TetEl)*cos(Phi),dyEl= sin(TetEl)*sin(Phi),dzEl=cos(TetEl);
  Float dxPo=-sin(TetPo)*cos(Phi),dyPo=-sin(TetPo)*sin(Phi),dzPo=cos(TetPo);
   
  //
  // kinematic of the created pair
  //
  // the electron and positron are assumed to have a symetric
  // angular distribution with respect to the Z axis along the parent photon.

  Float temp = 0.0;
  Float ElectKineEnergy = max(temp,ElectTotEnergy - kElectronMass);

  GPUThreeVector ElectDirection;
  ThreeVector_Set(ElectDirection,dxEl,dyEl,dzEl);
  ThreeVector_Rotate(ElectDirection,GammaDirection);
  ThreeVector_Extend(ElectDirection,sqrt(ElectTotEnergy*ElectTotEnergy-kElectronMass*kElectronMass));  

  // create G4DynamicParticle object for the particle1  
  // G4DynamicParticle* aParticle1= new G4DynamicParticle(
  //        theElectron,ElectDirection,ElectKineEnergy);
  GPUTrack *aParticle1 = &tracks[child_index];
  aParticle1->particle_id = 11;
  aParticle1->particle_index = electron_index;
  aParticle1->charge = -1;
  FourVector_Set(aParticle1->momentum,ElectDirection,ElectTotEnergy);
  
  // the e+ is always created (even with Ekine=0) for further annihilation.

  temp = 0.0;
  Float PositKineEnergy = max(temp,PositTotEnergy - kElectronMass);

  GPUThreeVector PositDirection;
  ThreeVector_Set(PositDirection,dxPo,dyPo,dzPo);
  ThreeVector_Rotate(PositDirection,GammaDirection);  
  ThreeVector_Extend(PositDirection,sqrt(PositTotEnergy*PositTotEnergy-kElectronMass*kElectronMass)); 

  // create G4DynamicParticle object for the particle2 
  // G4DynamicParticle* aParticle2= new G4DynamicParticle(
  //                     thePositron,PositDirection,PositKineEnergy);
  GPUTrack *aParticle2 = &tracks[child_index - 1];
  aParticle2->particle_id = 11;
  aParticle2->particle_index = electron_index;
  aParticle2->charge = 1;
  FourVector_Set(aParticle2->momentum,PositDirection,PositTotEnergy);

  // Kill photon
  CUDA_SetEnergy(track->momentum,0.0,0.0,index);

  // Spawn children
  SpawnChild(*track,child_index,kElectronMass);
  SpawnChild(*track,child_index-1,kElectronMass);
}

} // End unnamed namespace

__device__
void CUDA_GEANT4PairProduction(
    GPUTrack* track,
    const ParticlePars* particle,
    const MaterialPars* material,
    const Float dl,
    curandState *rng_state, 
    const int index) {

  if (track->particle_id != 22) return;

  int child_index = CanHaveTwoChildren(index);
  // Must be able to spawn two children
  if (child_index == -1) {
    UpdateState(index,WAITING);
    return;
  }

  // Use radiation length for probability
  Float chance_to_interact = 1 - exp(-dl/material->radiation_length);
  if (curand_uniform(rng_state) > chance_to_interact) return;

  SampleSecondaries(track,particle,material,index,child_index,rng_state);

}

} // End namespace na63