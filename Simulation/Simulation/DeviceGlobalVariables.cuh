#ifndef NA63_DEVICEGLOBALVARIABLES_CUH
#define NA63_DEVICEGLOBALVARIABLES_CUH

namespace na63 {

extern __constant__ VolumePars *volumes;
extern __constant__ int n_volumes;
extern __device__ GPUTrack *tracks;
extern __constant__ int maximum_index;
extern __device__ int *keys;
extern __device__ Float secondary_threshold;

extern __constant__ int electron_index;
extern __constant__ int photon_index;

}

#endif /* NA63_DEVICEGLOBALVARIABLES_CUH */