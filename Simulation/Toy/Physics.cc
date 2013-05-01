#include "Physics.hh"
#include <cmath>
#include <iostream>

// K/A = 0.307075 MeV g^-1 cm^2
// Z_FE = 26
// m_e = 0.51098892 MeV
// I_FE = 286

void BetheIron(Track* track, Float dl) {
  Float beta_squared = pow(track->beta(),2);
  Float gamma = track->gamma();
  Float gamma_squared = pow(gamma,2);
  Float mass = track->m();
  // PDG p. 287
  Float T_max = 9.18507876e16 * beta_squared * gamma_squared /
          (1 + 1.02197784 * gamma / mass + 0.261109676 / pow(mass,2));
  // PDG p. 286
  Float dE = dl * 26 * 0.307075 *
             pow(track->q(),2) / beta_squared;
  dE *= 1/2 * log(1.12292518e12 * beta_squared * gamma_squared * T_max)
        - beta_squared;
  std::cout << track->momentum[3] << " - " << dE << " = " << track->momentum[3] - dE << std::endl;
  track->momentum[3] -= dE;
}