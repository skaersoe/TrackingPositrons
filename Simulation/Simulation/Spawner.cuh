#ifndef NA63_SIMULATION_SPAWNER_CUH
#define NA63_SIMULATION_SPAWNER_CUH

#include "Geometry/LibraryCUDA.cuh"
#include "Simulation/DeviceGlobalVariables.cuh"

namespace na63 {

__device__ inline
GPUTrack* CanHaveChild(const int index) {
  int child_index = maximum_index - 2*index;
  GPUTrack *t = &track_pool[child_index];
  if (t->state == STATE_FREE) return t;
  return nullptr;
}

__device__ inline
GPUTrack* CanHaveTwoChildren(const int index) {
  int child_index = maximum_index - 2*index;
  GPUTrack *t1 = &track_pool[child_index];
  GPUTrack *t2 = &track_pool[child_index-1];
  if (t1->state == STATE_FREE &&
      t2->state == STATE_FREE) return t1;
  child_index *= 2;
  return nullptr;
}

__device__ inline
void SpawnChild(const GPUTrack& parent, GPUTrack *dst, GPUTrack& track) {
  track.state = STATE_ALIVE;
  FourVector_Copy(track.position,parent.position);
  track.volume_index = parent.volume_index;
  track.initial_index = parent.initial_index;
  FourVector_Copy(track.vertex,track.position);
  *dst = track;
}

} // End namespace na63

#endif /* NA63_SIMULATION_SPAWNER_CUH */