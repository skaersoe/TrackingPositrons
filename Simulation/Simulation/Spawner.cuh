#ifndef NA63_SIMULATION_SPAWNER_CUH
#define NA63_SIMULATION_SPAWNER_CUH

#include "Geometry/LibraryCUDA.cuh"
#include "Simulation/DeviceGlobalVariables.cuh"

namespace na63 {

__device__ inline
int CanHaveChild(const int index) {
  int child_index = maximum_index - 2*index;
  if (child_index < 0) return -1;
  GPUTrack *t = &tracks[child_index];
  // if (keys[index] == TRACK_KEY_WAITING) {
  //   printf("Waiting key %i looking at index %i, ",index,child_index);
  // }
  if (t->state == FREE) {
    // printf("Index %i can spawn at %i\n",index,child_index);
    // printf(" and it's available!\n");
    return child_index;
  }
  // printf(" but it's unavailable with state %i\n",t->state);
  return -1;
}

__device__ inline
int CanHaveTwoChildren(const int index) {
  int child_index = maximum_index - 2*index;
  if (child_index < 0) return -1;
  GPUTrack *t1 = &tracks[child_index];
  GPUTrack *t2 = &tracks[child_index-1];
  // if (keys[index] == TRACK_KEY_WAITING) {
  //   printf("Waiting key %i looking at index %i and %i, ",index,child_index,child_index-1);
  // }
  if (t1->state == FREE &&
      t2->state == FREE) {
    // printf(" and it's available!\n");
    // printf("Index %i can spawn at %i and %i\n",index,child_index,child_index-1);
    return child_index;
  }
  // printf(" but it's unavailable with state %i and %i\n",t1->state,t2->state);
  return -1;
}

__device__ inline
void SpawnChild(const GPUTrack& parent, const int index,
    const Float mass = 0.0) {
  GPUTrack *track = &tracks[index];
  // if (track->momentum[3] != track->momentum[3]) {
  //   printf("SpawnChild received NaN momentum track\n");
  //   return;
  // }
  track->volume_index = parent.volume_index;
  track->initial_index = parent.initial_index;
  FourVector_Copy(track->position,parent.position);
  FourVector_Copy(track->vertex,parent.position);
  if (track->momentum[3] <= mass) {
    UpdateState(index,DEAD);
  } else {
    UpdateState(index,ALIVE);
  }
  // printf("Spawning at %i with id %i, momentum %g, time %g\n",index,track->particle_id,track->momentum[3],track->position[3]);
}

} // End namespace na63

#endif /* NA63_SIMULATION_SPAWNER_CUH */