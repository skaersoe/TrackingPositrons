#ifndef NA63_CONSTANTS_H
#define NA63_CONSTANTS_H

#include <cmath>

#include "Geometry/Library.hh"

namespace na63 {

const Float ns = 1.0;
const Float  s = 1e9 * ns;
const Float ms = 1e-3 * s;

const Float cm = 1.0;
const Float m = 100.0 * cm;
const Float mm = cm / 10;

const Float e_SI = 1.60217733e-19;

const Float MeV = 1.0;
const Float keV = 1e3 * MeV;
const Float eV = 1e6 * MeV;
const Float joule = eV / e_SI;

const Float eplus = 1.0;
const Float e_squared = eplus * eplus;

const Float Megavolt = MeV / eplus;
const Float volt = 1e-6 * Megavolt;

const Float coulomb = eplus / e_SI;
const Float ampere = coulomb / s;
const Float weber = volt * s;
const Float henry = weber / ampere;

const Float kPi = 3.14159265;
const Float kC = 2.99792458e8 * m/s;
const Float kCSquared = kC * kC;
const Float kHPlanck = 6.6260755e-34 * joule*s;
const Float kHBarPlanck = kHPlanck / (2 * kPi);
const Float kHBarC = kHBarPlanck * kC;
const Float kElectronMass = 5.10998910e-1 * MeV;
const Float kMuonMass = 1.056583715e2 * MeV;
const Float fac_fel = log(184.15);
const Float fac_finel = log(1194.0);
const Float kMu0 = 4 * kPi * 1.0e-7 * henry / m;
const Float kEpsilon0 = 1.0 / (kCSquared * kMu0);
const Float kElmCoupling = e_squared / (4 * kPi * kEpsilon0);
const Float kElectronRadius = kElmCoupling / kElectronMass;
const Float kElectronComptonLength = kHBarC / kElectronMass;
const Float kMigdalConstant = kElectronRadius * kElectronComptonLength
    * kElectronComptonLength * 4.0 * kPi;
const Float kFineStructure = kElmCoupling / kHBarC;
const Float kBremFactor = kFineStructure * kElectronRadius*kElectronRadius
    * 16.0/3.0;

} // End namespace na63

#endif /* NA63_CONSTANTS_H */