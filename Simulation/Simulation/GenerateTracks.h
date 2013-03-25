#ifndef NA63_SIMULATION_GENERATEPARTICLES_H
#define NA63_SIMULATION_GENERATEPARTICLES_H

// #include "HepMC/GenParticle.h"
#include <cuda.h>
#include <cuda_runtime.h>

#include "Simulation/Simulator.h"
#include "Simulation/Track.h"

namespace na63 {

  void GenerateTracks(Track *t, SimulatorPars args);

}

#endif