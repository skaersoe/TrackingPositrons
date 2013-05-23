#ifndef NA63_TYPES
#define NA63_TYPES

#include <vector>


namespace na63 {

  // A detector hit.
  // Pixel coordinate, coordinate and detector number.

  typedef struct {
    int x_pixel;
    int y_pixel;
    int z_pixel;
    double x;
    double y;
    double z;
    int det;
  } hit;

  // Particle state.

  typedef struct {
    double x;
    double y;
    double z;
    double dx;
    double dy;
    double p;
    int q;
  } pstate;

  // Initial track.

  typedef struct {
    std::vector hits;
    double dx;
    double dy;
    double p;
    int q; 
  } track;

}

#endif