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
    externalParticles = true;
  }
  Simulator::~Simulator() {
    deleteParticles();
  }

  simulator_args_t Simulator::getArgs() {
    return args;
  }

  void Simulator::setArgs(simulator_args_t args) {
    this->args = args;
  }
  
  void Simulator::deleteParticles() {
    if (!externalParticles) {
      delete particles;
      args.N = 0;
    }
  }

  void Simulator::setParticles(simple_particle_t *particles, const unsigned N) {
    deleteParticles();
    args.N = N;
    this->particles = particles;
  }

  void Simulator::generateParticles() {
    deleteParticles();
    particles = new simple_particle_t[args.N];
    Simulation::generateParticles(particles,args);
    externalParticles = false;
  }

  void Simulator::propagate() {
    Simulation::propagate(particles,args);
  }

}