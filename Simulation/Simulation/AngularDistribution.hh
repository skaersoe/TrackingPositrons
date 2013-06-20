#ifndef NA63_SIMULATION_ANGULARDISTRIBUTION_H
#define NA63_SIMULATION_ANGULARDISTRIBUTION_H

#include "Geometry/Library.hh"

namespace na63 {

/** Generates random angular distributions for bremsstrahlung */

ThreeVector ModifiedTsai_SampleDirection(const Track& track);

} // End namespace na63

#endif /* NA63_SIMULATION_ANGULARDISTRIBUTION_H */