#include "Geometry/SimpleBox.hh"

namespace na63 {

bool SimpleBox::Inside(const FourVector& position) const {
  return position[0] <= center[0] + 0.5*size[0] &&
         position[0] >= center[0] - 0.5*size[0] &&
         position[1] <= center[1] + 0.5*size[1] &&
         position[1] >= center[1] - 0.5*size[1] &&
         position[2] <= center[2] + 0.5*size[2] &&
         position[2] >= center[2] - 0.5*size[2];
}

void SimpleBox::SetSpecificParameters(void *parameters) {
  SimpleBoxPars* simple_box = (SimpleBoxPars*)parameters;
  center.GPU(simple_box->center);
  size.GPU(simple_box->size);
}


} // End namespace na63