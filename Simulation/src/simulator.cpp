#include <iostream>

#include "Simulation/Simulator.h"
#include "Simulation/GenerateParticles.h"
#include "Simulation/Propagate.h"

namespace Simulation {

  Simulator::Simulator(void) : Simulator::Simulator(args) {
    simulator_args_t args = {
      .device = GPU,
      .debug = false
    };
  }
  Simulator::Simulator(simulator_args_t args) {
    setArgs(args);
    particles = NULL;
  }
  Simulator::~Simulator() {
    delete particles;
  }

  simulator_args_t Simulator::getArgs() {
    return args;
  }

  void Simulator::setArgs(simulator_args_t args) {
    this->args = args;
  }

  void Simulator::generateParticles() {
    delete particles;
    particles = new simple_particle_t[args.nParticles];
    Simulation::generateParticles(particles,args);
  }

  void Simulator::propagate() {
    Simulation::propagate(particles,args);
  }

}