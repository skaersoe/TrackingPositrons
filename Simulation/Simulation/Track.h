#ifndef NA63_SIMULATION_TRACK_H
#define NA63_SIMULATION_TRACK_H

namespace na63 {

  typedef struct {
    int particle_index; // Index into Particle array, managed dynamically
    float r[3];
    float p[3];
    // Align to 32 bytes
    char padding[4];
  } Track;

}

#endif /* NA63_SIMULATION_TRACK_H */