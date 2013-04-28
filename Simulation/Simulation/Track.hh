#ifndef NA63_SIMULATION_TRACK_H
#define NA63_SIMULATION_TRACK_H

namespace na63 {

  /**
   * Track data struct, should always be aligned to 32 bytes to take advantage
   * of memory bursts.
   */
  typedef struct {
    // int particle_index;
    int particle_id; // According to Monte Carlo Particle Numbering Scheme
                     // Will be changed to internal index at propagation time
    float r[3];
    float p[3];
    // Align to 32 bytes
    char padding[4];
  } Track;

}

#endif /* NA63_SIMULATION_TRACK_H */