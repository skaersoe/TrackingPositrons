#include "Volume.hh"

#ifndef GEOMETRY_BOX_CUH
#define GEOMETRY_BOX_CUH

namespace na63 {

  typedef struct {
  ThreeVector center;
  ThreeVector x_vector;
  ThreeVector y_vector;
  ThreeVector z_vector;
  // Pad to volume parameter size
  char padding[(int)(
    VOLUME_PARAMETER_SIZE
    - 4*sizeof(ThreeVector)
  )];
} BoxPars;

  __device__ inline
  bool Box_InsideKernel(ThreeVector point, void* pars);

  __device__ inline
  ThreeVector VectorAdd(ThreeVector a, ThreeVector b);

  __device__ inline
  ThreeVector VectorSubtract(ThreeVector a, ThreeVector b);

  __device__ inline
  float VectorDot(ThreeVector a, ThreeVector b);

  __device__ inline
  ThreeVector VectorNegate(ThreeVector a);

}

#endif /* GEOMETRY_BOX_CUH */