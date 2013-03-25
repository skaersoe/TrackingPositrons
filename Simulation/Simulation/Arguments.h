#ifndef NA63_SIMULATION_LAUNCHARGUMENTS_H
#define NA63_SIMULATION_LAUNCHARGUMENTS_H

namespace na63 {

  typedef enum {
    CPU,
    GPU
  } SimulatorDevice;

  typedef struct {
    SimulatorDevice device;
    bool debug;
    bool render;
    unsigned N;
    float dt;
    int steps;
  } SimulatorPars;

}

#endif