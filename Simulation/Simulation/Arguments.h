#ifndef LAUNCHARGUMENTS_H
#define LAUNCHARGUMENTS_H

namespace Simulation {

  typedef enum {
    CPU,
    GPU
  } simulator_device_t;

  typedef struct {
    simulator_device_t device;
    bool debug;
    unsigned nParticles;
  } simulator_args_t;

}

#endif