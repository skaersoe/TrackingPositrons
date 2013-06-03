#ifndef NA63_DEVICEGLOBALVARIABLES_CUH
#define NA63_DEVICEGLOBALVARIABLES_CUH

namespace na63 {

extern __constant__ VolumePars *volumes;
extern __constant__ int n_volumes;
extern __constant__ GPUTrack *track_pool;
extern __constant__ int maximum_index;

extern __constant__ int electron_index;
extern __constant__ int photon_index;

}

#endif /* NA63_DEVICEGLOBALVARIABLES_CUH */