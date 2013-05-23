#ifndef NA63_TRACKSEEDER
#define NA63_TRACKSEEDER

#include <vector>
#include "FileHandler.hh"

namespace na63 {

  // Contains the initial parameters of a seeded track.
  // Also contains a vector with the hits of the track. 

  typedef struct {
    std::vector hits;
    double x;
    double y;
    double dx;
    double dy;
    double p; 
  } track;

  // Uses the chi-square method of seeding tracks,
  // by taking each hit and matching it to a straight
  // line in the plane normal to the B-field direction.

 std::vector<seed> Seed(std::vector<Hit> hits);

}

#endif