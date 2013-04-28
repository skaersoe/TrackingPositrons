#include "../Geometry/Box.cuh"

#define SIGN_MASK 1 << sizeof(float)*8 - 1

namespace na63 {

  __device__ 
  bool Box_InsideKernel(ThreeVector point, void* pars) {

    BoxPars* boxpars = (BoxPars*)pars;

    ThreeVector x_top;
    ThreeVector x_bottom;
    ThreeVector x_negative;
    float x_top_dot;
    float x_bottom_dot;
    ThreeVector x_top_distance;
    ThreeVector x_bottom_distance;

    ThreeVector y_top;
    ThreeVector y_bottom;
    ThreeVector y_negative;
    float y_top_dot;
    float y_bottom_dot;
    ThreeVector y_top_distance;
    ThreeVector y_bottom_distance;

    ThreeVector z_top;
    ThreeVector z_bottom;
    ThreeVector z_negative;
    float z_top_dot;
    float z_bottom_dot;
    ThreeVector z_top_distance;
    ThreeVector z_bottom_distance;

    x_top = VectorAdd(boxpars -> center, boxpars -> x_vector);
    x_bottom = VectorSubtract(boxpars -> center, boxpars -> x_vector);
    x_top_distance = VectorSubtract(point, x_top);
    x_top_dot = VectorDot(boxpars -> x_vector, x_top_distance);
    x_bottom_distance = VectorSubtract(point, x_bottom);
    x_negative = VectorNegate(boxpars -> x_vector);
    x_bottom_dot = VectorDot(x_negative, x_bottom_distance);

    y_top = VectorAdd(boxpars -> center, boxpars -> y_vector);
    y_bottom = VectorSubtract(boxpars -> center, boxpars -> y_vector);
    y_top_distance = VectorSubtract(point, y_top);
    y_top_dot = VectorDot(boxpars -> y_vector, y_top_distance);
    y_bottom_distance = VectorSubtract(point, y_bottom);
    y_negative = VectorNegate(boxpars -> y_vector);
    y_bottom_dot = VectorDot(y_negative, y_bottom_distance);

    z_top = VectorAdd(boxpars -> center, boxpars -> z_vector);
    z_bottom = VectorSubtract(boxpars -> center, boxpars -> z_vector);
    z_top_distance = VectorSubtract(point, z_top);
    z_top_dot = VectorDot(boxpars -> z_vector, z_top_distance);
    z_bottom_distance = VectorSubtract(point, z_bottom);
    z_negative = VectorNegate(boxpars -> z_vector);
    z_bottom_dot = VectorDot(z_negative, z_bottom_distance);

    if ((x_top_dot & sign_mask) == (x_bottom_dot & sign_mask))
      return false;

    if ((y_top_dot & sign_mask) == (y_bottom_dot & sign_mask))
      return false;

    if ((z_top_dot & sign_mask) == (z_bottom_dot & sign_mask))
      return false;

    return true;

  }

  __device__ inline
  ThreeVector VectorAdd(ThreeVector a, ThreeVector b) {
    ThreeVector c;

    c.x = a.x + b.x;
    c.y = a.y + b.y;
    c.z = a.z + b.z;

    return c;
  }

  __device__ inline
  ThreeVector VectorSubtract(ThreeVector a, ThreeVector b) {
    ThreeVector c;

    c.x = a.x - b.x;
    c.y = a.y - b.y;
    c.z = a.z - b.z;

    return c;
  }

  __device__ inline
  float VectorDot(ThreeVector a, ThreeVector b) {
    float result;

    result = a.x*b.x + a.y*b.y + a.z*b.z;

    return result;
  }

  __device__ inline
  ThreeVector VectorNegate(ThreeVector a) {
    ThreeVector result;

    result.x = -a.x;
    result.y = -a.y;
    result.z = -a.z;

    return result;
  }

}