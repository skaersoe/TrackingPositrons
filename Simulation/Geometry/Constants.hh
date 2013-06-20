#ifndef NA63_CONSTANTS_H
#define NA63_CONSTANTS_H

#include <cmath>

#include "Geometry/Library.hh"

namespace na63 {

const Float ns = 1.0;
const Float  s = 1e9 * ns;
const Float ms = 1e-3 * s;

const Float cm = 1.0;
const Float m = 100 * cm;
const Float mm = 0.1 * cm;

const Float e_SI = 1.60217733e-19;

const Float MeV = 1.0;
const Float GeV = 1e3 * MeV;
const Float keV = 1e-3 * MeV;
const Float eV = 1e-6 * MeV;
const Float joule = eV / e_SI;

const Float eplus = 1.0;
const Float e_squared = eplus * eplus;

const Float Megavolt = MeV / eplus;
const Float volt = 1e-6 * Megavolt;

const Float coulomb = eplus / e_SI;
const Float ampere = coulomb / s;
const Float weber = volt * s;
const Float henry = weber / ampere;

const Float mol = 1;

const Float kPi = 3.14159265;
const Float kC = 2.99792458e8 * m/s;
const Float kCSquared = kC * kC;
const Float kHPlanck = 6.6260755e-34 * joule*s;
const Float kHBarPlanck = kHPlanck / (2 * kPi);
const Float kHBarC = kHBarPlanck * kC;
const Float kElectronMass = 0.510998910 * MeV;
const Float kMu0 = 4 * kPi * 1.0e-7 * henry / m;
const Float kEpsilon0 = 1.0 / (kCSquared * kMu0);
const Float kElmCoupling = e_squared / (4 * kPi * kEpsilon0);
const Float kElectronRadius = kElmCoupling / kElectronMass;
const Float kElectronRadiusSquared = kElectronRadius*kElectronRadius;
const Float kElectronComptonLength = kHBarC / kElectronMass;
const Float kMigdalConstant = kElectronRadius * kElectronComptonLength
    * kElectronComptonLength * 4.0 * kPi;
const Float kFineStructure = kElmCoupling / kHBarC;
const Float kBremFactor = kFineStructure * kElectronRadius*kElectronRadius
    * 16.0/3.0;
const Float kAvogadro = 6.022e23 * 1/mol; // Mol^-1

const Float kMuonMass = 1.05658372e2 * MeV; // MeV
const Float kMuonCharge = -1.0 * eplus;
const Float kProtonCharge = 1.0 * eplus;
const Float kProtonMass = 938.272046 * MeV;
// http://www.physics.nist.gov/cgi-bin/Star/compos.pl?refer=ap&matno=026
const Float kIronMeanExcitationPotential = 286.0e-6 * MeV;
const Float kIronAtomicNumber = 26.0;
const Float kIronDensity = 7.874; // g/cm^-3
const Float kIronAtomicWeight = 55.84;
const Float kIronRadiationLength = 13.84/kIronDensity;
const Float kCopperAtomicNumber = 29.0;
const Float kCopperDensity = 8.96;
const Float kCopperAtomicWeight = 63.546;
const Float kCopperMeanExcitationPotential = 322.0e-6 * MeV;
const Float kCopperRadiationLength = 12.86/kCopperDensity;
const Float kLeadAtomicNumber = 82.0;
const Float kLeadDensity = 11.35;
const Float kLeadRadiationLength = 0.5612;
const Float kLeadAtomicWeight = 207.2;
const Float kLeadMeanExcitationPotential = 823.0e-6;
const Float kGraphiteAtomicNumber = 6.0;
const Float kGraphiteDensity = 2.21;
const Float kGraphiteRadiationLength = 19.32;
const Float kGraphiteAtomicWeight = 12.0107;
const Float kGraphiteMeanExcitationPotential = 78.0e-6;
const Float kLPMConstant = kFineStructure*kElectronMass*kElectronMass/(4.0*kPi*kHBarC)*0.5;
#ifndef __CUDACC_
const Float kLogTwo = log(2.0);
const Float fac_fel = log(184.15);
const Float fac_finel = log(1194.0);
#else
const Float kLogTwo = logf(2.0);
const Float fac_fel = logf(184.15);
const Float fac_finel = logf(1194.0);
#endif

} // End namespace na63

#endif /* NA63_CONSTANTS_H */