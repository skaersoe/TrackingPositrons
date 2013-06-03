#include <cmath>
#include <algorithm>

#include "Simulation/Bremsstrahlung.hh"
#include "Simulation/AngularDistribution.hh"
#include "Geometry/Library.hh"
#include "Geometry/Constants.hh"

namespace na63 {

Bremsstrahlung::Bremsstrahlung(const Particle* p) :
    particle(nullptr),
    is_electron(true) {

  min_threshold = 0.1 * keV;
  low_kinetic_energy = 10.0 * MeV;
  low_limit = 10.0 * MeV;

  particle_mass = kinetic_energy = total_energy = Z = z13 = z23 = lnZ = fel
      = finel = density_factor = f_max = coulomb_correction = 0;

  if (p) SetParticle(p);

}
    
/** Used in SampleSecondaries() */
void Bremsstrahlung::SetMaterial(const Material* material_new) {
  
  // std::cout<<"SetCurrentElement Z="<<Z<<std::endl;

  if (material != material_new) {

    material = material_new;
    Z = material->atomic_number();

    int iz = int(Z);
    // z13 = nist->GetZ13(iz);
    z13 = pow(iz, 0.33333333333333);
    z23 = z13*z13;
    
    lnZ = log(iz);
    // lnZ = nist->GetLOGZ(iz);

    fel = fac_fel - lnZ/3.0;
    finel = fac_finel - 2.0*lnZ/3.0;

    coulomb_correction = material->coulomb_correction();

    /*
    printf("------------------------\n");
    printf("fCoulomb :%f\n", fCoulomb);
    printf("------------------------\n");
    */
    
    f_max = fel-coulomb_correction + finel/Z + (1.0 + 1.0/Z)/12.0;

  }

} // End Bremsstrahlung::SetMaterial

/** Used in constructor */
void Bremsstrahlung::SetParticle(const Particle* particle_new) {

  if (particle == particle_new) return;

  particle = particle_new;
  particle_mass = particle->mass();

  if (particle->id() == 11) {
    is_electron = true;
  } else {
    is_electron = false;
  }

} // End Bremsstrahlung::SetParticle

/** Used in SampleSecondaries() */
void Bremsstrahlung::SetupForMaterial(const Material* mat,
    const Float kin_energy) {

  density_factor = mat->electron_density() * kMigdalConstant;

  // Calculate threshold for density effect
  kinetic_energy = kin_energy;
  total_energy = kinetic_energy + particle_mass;
  density_correction = density_factor * total_energy * total_energy; 

}

/** Used in SampleSecondaries(),
            ComputeCrossSectionPerAtom() */
Float Bremsstrahlung::ComputeDXSectionPerAtom(const Float gamma_energy) const {

  if (gamma_energy < 0.0) return 0.0;

  Float y = gamma_energy / total_energy;

  Float main = 0.0;
  //secondTerm=0.;

  // ** form factors complete screening case **
  // only valid for high energies (and if LPM suppression does not play a role)
  main = (3.0 / 4.0 * y*y - y + 10.0) * ((fel - coulomb_correction) + finel/Z);
  // secondTerm = (1.-y)/12.*(1.+1./currentZ);

  /*
  std::cout<<" F1(0) "<<ScreenFunction1(0.) <<std::endl;
  std::cout<<" F1(0) "<<ScreenFunction2(0.) <<std::endl;
  std::cout<<"Ekin = "<<kinEnergy<<std::endl;
  std::cout<<"Z = "<<currentZ<<std::endl;
  std::cout<<"main  = "<<main<<std::endl;
  std::cout<<" y = "<<y<<std::endl;
  std::cout<<" Fel-fCoulomb "<< (Fel-fCoulomb) <<std::endl;
  */

  Float main2 = ComputeParametrizedDXSectionPerAtom(kinetic_energy,gamma_energy,Z);

  /*
  std::cout<<"main2 = "<<main2<<std::endl;
  std::cout<<"main2tot = "<<main2 * ( (Fel-fCoulomb) + Finel/currentZ ) /  (Fel-fCoulomb);
  */

  Float cross = main2; //main+secondTerm;

  return cross;

} // End Bremsstrahlung::ComputeDXSectionPerAtom()

/** Used in ComputeCrossSectionPerAtom() */
Float Bremsstrahlung::ComputeXSectionPerAtom(const Float cut) const {

  Float cross = 0.0;

  // number of intervals and integration step 
  Float vcut = log(cut/total_energy);
  Float vmax = log(kinetic_energy/total_energy);
  int n = (int)(0.45*(vmax - vcut)) + 4;
  Float delta = (vmax - vcut)/Float(n);

  Float e0 = vcut;
  Float xs; 

  // integration
  for (int l=0;l<n;l++) {

    for (int i=0;i<8;i++) {

      Float eg = exp(e0 + xgi[i] * delta) * total_energy;
      xs = ComputeDXSectionPerAtom(eg);

      cross += wgi[i] * xs / (1.0 + density_correction / (eg*eg));

    }

    e0 += delta;
  }

  cross *= delta;

  return cross;

} // End Bremsstrahlung::ComputeXSectionPerAtom()

/** Used in ComputeCrossSectionPerAtom() */
Float Bremsstrahlung::ComputeCrossSectionPerAtom(const Float cut_energy, 
    const Float max_energy) const {

  if (kinetic_energy < low_kinetic_energy) return 0.0;

  Float cut  = std::min(cut_energy, kinetic_energy);
  Float tmax = std::min(max_energy, kinetic_energy);

  if (cut >= tmax) return 0.0;

  Float cross = ComputeXSectionPerAtom(cut);

  // Allow partial integration
  if (tmax < kinetic_energy) {
    cross -= ComputeXSectionPerAtom(tmax);
  }

  cross *= Z*Z*kBremFactor;

  return cross;

} // End Bremsstrahlung::ComputeCrossSectionPerAtom()

void Bremsstrahlung::SampleSecondaries(
    Simulator *simulator, 
    // MDJ: replaces MaterialCutsCouple kinda wrong, be aware...
    const Material* couple,
    Track& track,
    Float cut_energy,
    Float max_energy) {

  Float kin_energy = track.kinetic_energy();
  if (kin_energy < low_kinetic_energy) return;

  Float cut  = std::min(cut_energy, kin_energy);
  Float emax = std::min(max_energy, kin_energy);

  /*
  printf("eKin part: %f\n", kineticEnergy);
  printf("lowKinThreshold %f\n", lowKinEnergy);
  printf("cut %f\n",cut);
  printf("emax %f\n", emax);
  */

  if (cut >= emax) return;

  SetupForMaterial(couple,kin_energy);

  // in VEmModel.cc get element based on cross section
  //const Element* elm = SelectRandomAtom(couple,particle,kineticEnergy,cut,emax);
  SetMaterial(material);

  Float xmin = log(cut*cut + density_correction);
  Float xmax = log(emax*emax  + density_correction);
  Float gamma_energy, f, x; 

  do {
    x = exp(xmin + gRandom->Uniform()*(xmax - xmin)) - density_correction;
    if (x < 0.0) x = 0.0;
    gamma_energy = sqrt(x);
    f = ComputeDXSectionPerAtom(gamma_energy);

    if (f > f_max) {
      std::cout << "### eBremParametrizedModel Warning: Majoranta exceeded! "
                << f << " > " << f_max
                << " Egamma(MeV)= " << gamma_energy
                << " E(mEV)= " << kinetic_energy
                << std::endl;
    }
    // printf("f %f\n", f);

  } while (f < f_max * gRandom->Uniform());

  // Angles of the emitted gamma. ( Z - axis along the parent particle)
  // Use general interface
  ThreeVector gamma_momentum = ModifiedTsai_SampleDirection(track);
  gamma_momentum.Extend(sqrt(gamma_energy));
  Track gamma(22,0,track.position,FourVector(gamma_momentum[0],gamma_momentum[1],gamma_momentum[2],gamma_energy));

  simulator->AddTrack(gamma);

  Float total_momentum = sqrt(kinetic_energy*(total_energy + kElectronMass));
  ThreeVector direction = track.momentum.Normalized();
  direction.Extend(total_momentum);
  direction -= gamma_momentum;
  direction.Normalize();

  // Energy and momentum of primary
  Float final_energy = kinetic_energy - gamma_energy;
  direction.Extend(final_energy*final_energy - kElectronMass*kElectronMass);
  FourVector momentum(direction[0],direction[1],direction[2],final_energy);

  if (gamma_energy > secondary_threshold) {

    // Stop tracking and create new secondary instead of primary
    
    /*  
    fParticleChange->ProposeTrackStatus(fStopAndKill);
    fParticleChange->SetProposedKineticEnergy(0.0);
    */

    track.Kill();

    Track electron(
      11,
      track.charge(),
      track.position,
      momentum
    );
    
    simulator->AddTrack(electron);

  } else {
    
    // Just update momentum and energy...

    track.momentum = momentum;

  }

}

Float ScreenFunction1(Float screen_variable) {

  // compute the value of the screening function 3*PHI1 - PHI2

  Float screen_value;

  if (screen_variable > 1.0) {
    screen_value = 42.24 - 8.368 * log(screen_variable + 0.952);
  } else {
    screen_value = 42.392 - screen_variable * (7.796 - 1.961 * screen_variable);
  }

  return screen_value;
} 

Float ScreenFunction2(Float screen_variable) {

  // compute the value of the screening function 1.5*PHI1 - 0.5*PHI2

  Float screen_value;

  if (screen_variable > 1.0) {
    screen_value = 42.24 - 8.368 * log(screen_variable + 0.952);
  } else {
    screen_value = 41.734 - screen_variable * (6.484 - 1.250 * screen_variable);
  }

  return screen_value;
}

/** Parametrized cross section.
    Used in SampleSecondaries() */
Float ComputeParametrizedDXSectionPerAtom(Float kinetic_energy, Float gamma_energy, Float Z) {

  static const double
     ah10 = 4.67733E+00, ah11 =-6.19012E-01, ah12 = 2.02225E-02,
     ah20 =-7.34101E+00, ah21 = 1.00462E+00, ah22 =-3.20985E-02,
     ah30 = 2.93119E+00, ah31 =-4.03761E-01, ah32 = 1.25153E-02;

  static const double
     bh10 = 4.23071E+00, bh11 =-6.10995E-01, bh12 = 1.95531E-02,
     bh20 =-7.12527E+00, bh21 = 9.69160E-01, bh22 =-2.74255E-02,
     bh30 = 2.69925E+00, bh31 =-3.63283E-01, bh32 = 9.55316E-03;

  static const double
     al00 =-2.05398E+00, al01 = 2.38815E-02, al02 = 5.25483E-04,
     al10 =-7.69748E-02, al11 =-6.91499E-02, al12 = 2.22453E-03,
     al20 = 4.06463E-02, al21 =-1.01281E-02, al22 = 3.40919E-04;

  static const double
     bl00 = 1.04133E+00, bl01 =-9.43291E-03, bl02 =-4.54758E-04,
     bl10 = 1.19253E-01, bl11 = 4.07467E-02, bl12 =-1.30718E-03,
     bl20 =-1.59391E-02, bl21 = 7.27752E-03, bl22 =-1.94405E-04;

  static const double t_low = 1.0;

  Float lnZ = log(Z); // 3.*(anElement->GetIonisation()->GetlogZ3());
  Float FZ = lnZ* (4.- 0.55*lnZ);
  Float ZZ = pow (Z*(Z+1.),1./3.); // anElement->GetIonisation()->GetZZ3();
  Float Z3 = pow (Z,1./3.); // (anElement->GetIonisation()->GetZ3())

  Float total_energy = kinetic_energy + kElectronMass;

  // Float x, epsil, greject, migdal, grejmax, q;
  Float epsil, greject;
  Float U  = log(kinetic_energy/kElectronMass);
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
    Float F1 = std::max(ScreenFunction1(screenvar) - FZ,Float(0.0));
    Float F2 = std::max(ScreenFunction2(screenvar) - FZ,Float(0.0));


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

} // End ComputeParametrizedDXSectionPerAtom()

} // End namespace na63