#ifndef SIMULATOR_H
#define SIMULATOR_H

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

      void generateParticles();
      void propagate();

    private:
      simulator_args_t args;
      simple_particle_t *particles;

  };

}

#endif