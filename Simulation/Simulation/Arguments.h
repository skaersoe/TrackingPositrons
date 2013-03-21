#ifndef NA63_SIMULATION_LAUNCHARGUMENTS_H
#define NA63_SIMULATION_LAUNCHARGUMENTS_H

namespace NA63 {

  typedef enum {
    CPU,
    GPU
  } simulator_device_t;

  typedef struct {
    simulator_device_t device;
    bool debug;
    bool render;
    unsigned N;
    float dt;
    int steps;
  } simulator_args_t;

}

#endif