#ifndef NA63_SIMULATION_BREMSSTRAHLUNG_H
#define NA63_SIMULATION_BREMSSTRAHLUNG_H

#include "Geometry/Library.hh"
#include "Simulation/Track.hh"
#include "Simulation/Simulator.hh"

namespace na63 {

void InitializeBremsstrahlung();

void Bremsstrahlung(Track& track, const Material& material,
    const Float dl);

} // End namespace na63

#endif /* NA63_SIMULATION_BRESSTRAHLUNG_H */