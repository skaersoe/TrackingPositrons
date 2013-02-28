#include <iostream>
#include "Simulation/Step.h"

void propagate(simple_particle_t *p, const int N) {
  std::cout << "Particle 5 was at ("
    << p[5].r[0] << ","
    << p[5].r[1] << ","
    << p[5].r[2] << ")...";
  for (int i=0;i<N;i++)
    for (int j=0;j<1000;j++)
      timestep(&p[i],0.1);
  std::cout << "...and is now at ("
    << p[5].r[0] << ","
    << p[5].r[1] << ","
    << p[5].r[2] << ").";
}