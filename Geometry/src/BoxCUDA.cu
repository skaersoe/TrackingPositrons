#include "Geometry/Box.hh"
#include "Geometry/BoxCUDA.cuh"
#include "Geometry/LibraryCUDA.cuh"

#define SIGN_MASK 1 << sizeof(Float)*8 - 1

namespace na63 {

__device__ 
bool Box_Inside(const GPUFourVector& point, const void* pars) {

  BoxPars* boxpars = (BoxPars*)pars;

  GPUThreeVector x_top;
  GPUThreeVector x_bottom;
  GPUThreeVector x_negative;
  Float x_top_dot;
  Float x_bottom_dot;
  GPUThreeVector x_top_distance;
  GPUThreeVector x_bottom_distance;

  GPUThreeVector y_top;
  GPUThreeVector y_bottom;
  GPUThreeVector y_negative;
  Float y_top_dot;
  Float y_bottom_dot;
  GPUThreeVector y_top_distance;
  GPUThreeVector y_bottom_distance;

  GPUThreeVector z_top;
  GPUThreeVector z_bottom;
  GPUThreeVector z_negative;
  Float z_top_dot;
  Float z_bottom_dot;
  GPUThreeVector z_top_distance;
  GPUThreeVector z_bottom_distance;

  ThreeVector_Add(boxpars->center,boxpars->x_vector,x_top);
  ThreeVector_Subtract(boxpars->center,boxpars->x_vector,x_bottom);
  ThreeVector_Subtract(point,x_top,x_top_distance);
  x_top_dot = ThreeVector_Dot(boxpars->x_vector,x_top_distance);
  ThreeVector_Subtract(point,x_bottom,x_bottom_distance);
  ThreeVector_Negate(boxpars->x_vector,x_negative);
  x_bottom_dot = ThreeVector_Dot(x_negative, x_bottom_distance);

  ThreeVector_Add(boxpars->center,boxpars->y_vector,y_top);
  ThreeVector_Subtract(boxpars->center,boxpars->y_vector,y_bottom);
  ThreeVector_Subtract(point,y_top,y_top_distance);
  y_top_dot = ThreeVector_Dot(boxpars->y_vector,y_top_distance);
  ThreeVector_Subtract(point,y_bottom,y_bottom_distance);
  ThreeVector_Negate(boxpars->y_vector,y_negative);
  y_bottom_dot = ThreeVector_Dot(y_negative,y_bottom_distance);

  ThreeVector_Add(boxpars->center,boxpars->z_vector,z_top);
  ThreeVector_Subtract(boxpars->center,boxpars->z_vector,z_bottom);
  ThreeVector_Subtract(point,z_top,z_top_distance);
  z_top_dot = ThreeVector_Dot(boxpars->z_vector,z_top_distance);
  ThreeVector_Subtract(point,z_bottom,z_bottom_distance);
  ThreeVector_Negate(boxpars->z_vector,z_negative);
  z_bottom_dot = ThreeVector_Dot(z_negative,z_bottom_distance);

  if ((x_top_dot < 0) == (x_bottom_dot < 0))
    return false;

  if ((y_top_dot < 0) == (y_bottom_dot < 0))
    return false;

  if ((z_top_dot < 0) == (z_bottom_dot < 0))
    return false;

  return true;

} // End namespace na63

InsideFunction Box::inside_function_ = Box_Inside;

}