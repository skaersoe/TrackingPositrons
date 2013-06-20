#include "Geometry/LibraryCUDA.cuh"
#include "Geometry/Constants.hh"
#include "Simulation/ScreenFunctionsCUDA.cuh"
#include "Simulation/AngularDistributionCUDA.cuh"
#include "Simulation/Spawner.cuh"
#include "Simulation/TrackGPU.cuh"
#include <cfloat>

namespace na63 {

namespace {

__device__ __constant__ Float
   ah10 = 4.67733E+00, ah11 =-6.19012E-01, ah12 = 2.02225E-02,
   ah20 =-7.34101E+00, ah21 = 1.00462E+00, ah22 =-3.20985E-02,
   ah30 = 2.93119E+00, ah31 =-4.03761E-01, ah32 = 1.25153E-02;

__device__ __constant__ Float
   bh10 = 4.23071E+00, bh11 =-6.10995E-01, bh12 = 1.95531E-02,
   bh20 =-7.12527E+00, bh21 = 9.69160E-01, bh22 =-2.74255E-02,
   bh30 = 2.69925E+00, bh31 =-3.63283E-01, bh32 = 9.55316E-03;

__device__ __constant__ Float
   al00 =-2.05398E+00, al01 = 2.38815E-02, al02 = 5.25483E-04,
   al10 =-7.69748E-02, al11 =-6.91499E-02, al12 = 2.22453E-03,
   al20 = 4.06463E-02, al21 =-1.01281E-02, al22 = 3.40919E-04;

__device__ __constant__ Float
   bl00 = 1.04133E+00, bl01 =-9.43291E-03, bl02 =-4.54758E-04,
   bl10 = 1.19253E-01, bl11 = 4.07467E-02, bl12 =-1.30718E-03,
   bl20 =-1.59391E-02, bl21 = 7.27752E-03, bl22 =-1.94405E-04;

__device__ __constant__ Float t_low = 1.0;

__device__ __constant__ Float low_kinetic_energy = 10.0 * MeV;

__device__ __constant__ Float fac_fel_d = 5.21575064289;
__device__ __constant__ Float fac_finel_d = 7.08506429395;

__device__
Float ComputeParametrizedDXSectionPerAtom(Float kinetic_energy, Float gamma_energy, Float Z) {


  Float lnZ = logf(Z); // 3.*(anElement->GetIonisation()->GetlogZ3());
  Float FZ = lnZ* (4.- 0.55*lnZ);
  Float ZZ = pow(Float(Z*(Z+1.)),Float(1.0/3.0)); // anElement->GetIonisation()->GetZZ3();
  Float Z3 = pow(Z,Float(1.0/3.0)); // (anElement->GetIonisation()->GetZ3())

  Float total_energy = kinetic_energy + kElectronMass;

  // Float x, epsil, greject, migdal, grejmax, q;
  Float epsil, greject;
  Float U  = logf(kinetic_energy/kElectronMass);
  Float U2 = U*U;

  // Precalculated parameters
  Float ah, bh;

  if (kinetic_energy > t_low) {
       
    Float ah1 = ah10 + ZZ * (ah11 + ZZ * ah12);
    Float ah2 = ah20 + ZZ * (ah21 + ZZ * ah22);
    Float ah3 = ah30 + ZZ * (ah31 + ZZ * ah32);

    Float bh1 = bh10 + ZZ * (bh11 + ZZ * bh12);
    Float bh2 = bh20 + ZZ * (bh21 + ZZ * bh22);
    Float bh3 = bh30 + ZZ * (bh31 + ZZ * bh32);

    ah = 1.0  + (ah1*U2 + ah2*U + ah3) / (U2*U);
    bh = 0.75 + (bh1*U2 + bh2*U + bh3) / (U2*U);

    // Limit of the screening variable
    Float screenfac = 136.0*kElectronMass/(Z3*total_energy);

    // epsil = x*kinetic_energy/total_energy;
    epsil = gamma_energy/total_energy;
    Float screenvar = screenfac*epsil/(1.0-epsil);
    Float F1 = max(CUDA_ScreenFunction1(screenvar) - FZ,Float(0.0));
    Float F2 = max(CUDA_ScreenFunction2(screenvar) - FZ,Float(0.0));


    greject = (F1 - epsil* (ah*F1 - bh*epsil*F2))/8.0; //  1./(42.392 - FZ);

    /*
    std::cout << " yy = "<<epsil<<std::endl;
    std::cout << " F1/(...) "<<F1/(42.392 - FZ)<<std::endl;
    std::cout << " F2/(...) "<<F2/(42.392 - FZ)<<std::endl;
    std::cout << " (42.392 - FZ) " << (42.392 - FZ) <<std::endl;
    */

  } else { // kinetic_energy < t_low

    Float al0 = al00 + ZZ* (al01 + ZZ* al02);
    Float al1 = al10 + ZZ* (al11 + ZZ* al12);
    Float al2 = al20 + ZZ* (al21 + ZZ* al22);
 
    Float bl0 = bl00 + ZZ* (bl01 + ZZ* bl02);
    Float bl1 = bl10 + ZZ* (bl11 + ZZ* bl12);
    Float bl2 = bl20 + ZZ* (bl21 + ZZ* bl22);
 
    ah = al0 + al1*U + al2*U2;
    bh = bl0 + bl1*U + bl2*U2;

    Float x = gamma_energy/kinetic_energy;
    greject = (1.0 + x* (ah + bh*x));

    /*
    // Compute the maximum of the rejection function
    grejmax = max(1. + xmin* (ah + bh*xmin), 1.+ah+bh);
    Float xm = -ah/(2.*bh);
    if ( xmin < xm && xm < xmax) grejmax = max(grejmax, 1.+ xm* (ah + bh*xm));
    */
  }

 return greject;

}

__device__
void SampleSecondaries(
    GPUTrack* track,
    // MDJ: replaces MaterialCutsCouple kinda wrong, be aware...
    const ParticlePars* particle,
    const MaterialPars* couple,
    const int index,
    const int child_index,
    curandState *rng_state,
    Float cut_energy = 2*kElectronMass,
    Float max_energy = FLT_MAX) {

  Float kin_energy = track->momentum[3] - particle->mass;
  if (kin_energy < low_kinetic_energy) return;

  Float cut  = min(cut_energy, kin_energy);
  Float emax = min(max_energy, kin_energy);

  /*
  printf("eKin part: %f\n", kineticEnergy);
  printf("lowKinThreshold %f\n", lowKinEnergy);
  printf("cut %f\n",cut);
  printf("emax %f\n", emax);
  */

  if (cut >= emax) return;

  // CUDA_GEANT4Bremsstrahlung_SetupForMaterial(couple,kin_energy);

  Float density_factor = couple->electron_density * kMigdalConstant;

  // Calculate threshold for density effect
  Float kinetic_energy = kin_energy;
  Float total_energy = kinetic_energy + particle->mass;
  Float density_correction = density_factor * total_energy * total_energy; 

  // in VEmModel.cc get element based on cross section
  //const Element* elm = SelectRandomAtom(couple,particle,kineticEnergy,cut,emax);

  Float Z = couple->atomic_number;

  int iz = int(Z);
  // z13 = nist->GetZ13(iz);
  Float z13 = pow(iz, 0.33333333333333);
  // Float z23 = z13*z13;
  
  Float lnZ = logf(iz);
  // lnZ = nist->GetLOGZ(iz);

  Float fel = fac_fel_d - lnZ/3.0;
  Float finel = fac_finel_d - 2.0*lnZ/3.0;

  Float coulomb_correction = couple->coulomb_correction;

  /*
  printf("------------------------\n");
  printf("fCoulomb :%f\n", fCoulomb);
  printf("------------------------\n");
  */
  
  Float f_max = fel-coulomb_correction + finel/Z + (1.0 + 1.0/Z)/12.0;

  Float xmin = logf(cut*cut + density_correction);
  Float xmax = logf(emax*emax + density_correction);
  Float gamma_energy, f, x; 

  do {
    x = exp(xmin + curand_uniform(rng_state)*(xmax - xmin)) - density_correction;
    if (x < 0.0) x = 0.0;
    gamma_energy = sqrt(x);
    // f = CUDA_GEANT4Bremsstrahlung_ComputeDXSectionPerAtom(gamma_energy,total_energy,Z);

    if (gamma_energy < 0.0) {
      f = 0;
    } else {
      f = ComputeParametrizedDXSectionPerAtom(kinetic_energy,gamma_energy,Z);
    }

  } while (f < f_max * curand_uniform(rng_state));

  // Angles of the emitted gamma. ( Z - axis along the parent particle)
  // Use general interface
  GPUThreeVector gamma_momentum;
  CUDA_ModifiedTsai_SampleDirection(gamma_momentum,track,particle->mass,rng_state);
  GPUTrack *gamma = &tracks[child_index];
  gamma->particle_id = 22;
  gamma->particle_index = photon_index;
  gamma->charge = 0;
  FourVector_Set(gamma->momentum,gamma_momentum,gamma_energy);

  Float total_momentum = sqrt(kinetic_energy*(total_energy + kElectronMass));
  GPUThreeVector direction;
  ThreeVector_Normalized(direction,track->momentum);
  ThreeVector_Extend(direction,total_momentum);
  ThreeVector_Subtract(direction,gamma_momentum,direction);
  ThreeVector_Normalize(direction);

  // Energy and momentum of primary
  Float final_energy = total_energy - gamma_energy;
  ThreeVector_Extend(direction,sqrt(final_energy*final_energy - kElectronMass*kElectronMass));
  GPUFourVector momentum;
  FourVector_Set(momentum,direction,final_energy);

  if (gamma_energy > secondary_threshold) {

    // Stop tracking and create new secondary instead of primary
    
    /*  
    fParticleChange->ProposeTrackStatus(fStopAndKill);
    fParticleChange->SetProposedKineticEnergy(0.0);
    */

    GPUTrack *electron = &tracks[child_index-1];
    electron->particle_id = 11;
    electron->particle_index = electron_index;
    electron->charge = -1;
    FourVector_Copy(electron->momentum,momentum);

    SpawnChild(*track,child_index-1,kElectronMass);

    // Force track to kill itself...
    CUDA_SetEnergy(track->momentum,particle->mass,0.0,index);


  } else {
    
    // Just update momentum and energy...

    CUDA_SetMomentum(track->momentum,particle->mass,momentum,index);

  }

  SpawnChild(*track,child_index);

} // End SampleSecondaries()

} // End unnamed namespace

__device__
void CUDA_GEANT4Bremsstrahlung(
    GPUTrack* track,
    const ParticlePars* particle,
    const MaterialPars* material,
    const Float dl,
    curandState *rng_state, 
    const int index) {

  if (track->particle_id != 11) return;

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