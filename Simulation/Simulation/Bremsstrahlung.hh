#ifndef NA63_SIMULATION_BREMSSTRAHLUNG_H
#define NA63_SIMULATION_BREMSSTRAHLUNG_H

#include "Geometry/Library.hh"
#include "Simulation/Track.hh"
#include "Geometry/Material.hh"

namespace na63 {

void MyBremsstrahlung(Track& mother, const Material& material,
    const Float dl);

} // End namespace na63

#endif /* NA63_SIMULATION_BREMSSTRAHLUNG_H */