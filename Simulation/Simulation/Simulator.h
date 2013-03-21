#ifndef NA63_SIMULATION_SIMULATOR_H
#define NA63_SIMULATION_SIMULATOR_H

#include <thrust/host_vector.h>

#include "Simulation/Particle.h"
#include "Simulation/Arguments.h"

namespace Simulation {

  class Simulator {

    public:
      Simulator(void);
      Simulator(simulator_args_t args);

      ~Simulator();

      simulator_args_t getArgs();

      void setArgs(simulator_args_t args);
      void setParticles(simple_particle_t *particles, const unsigned N);

      void generateParticles();
      void propagate();

    private:
      simulator_args_t args;
      simple_particle_t *particles;
      bool externalParticles;

      void deleteParticles();

  };

}

#endif