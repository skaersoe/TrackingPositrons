#include "Geometry/SimpleBox.hh"
#include "Geometry/SimpleBoxCUDA.cuh"

namespace na63 {

__device__
bool SimpleBox_Inside(const GPUFourVector& point, const void* parameters) {
  SimpleBoxPars *simple_box = (SimpleBoxPars*)parameters;
  // printf("Checking if (%f,%f,%f) is inside (%f,%f,%f) with (%f,%f,%f)\n",
  //   point[0],
  //   point[1],
  //   point[2],
  //   simple_box->center[0],
  //   simple_box->center[1],
  //   simple_box->center[2],
  //   simple_box->size[0],
  //   simple_box->size[1],
  //   simple_box->size[2]
  // );
  return point[0] <= simple_box->center[0] + 0.5*simple_box->size[0] &&
         point[0] >= simple_box->center[0] - 0.5*simple_box->size[0] &&
         point[1] <= simple_box->center[1] + 0.5*simple_box->size[1] &&
         point[1] >= simple_box->center[1] - 0.5*simple_box->size[1] &&
         point[2] <= simple_box->center[2] + 0.5*simple_box->size[2] &&
         point[2] >= simple_box->center[2] - 0.5*simple_box->size[2];
}

InsideFunction SimpleBox::inside_function_ = SimpleBox_Inside;

} // End namespace na63