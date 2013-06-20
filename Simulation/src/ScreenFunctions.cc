#include "Geometry/Library.hh"

namespace na63 {

Float ScreenFunction1(Float screen_variable) {

  // compute the value of the screening function 3*PHI1 - PHI2

  Float screen_value;

  if (screen_variable > 1.0) {
    screen_value = 42.24 - 8.368 * log(screen_variable + 0.952);
  } else {
    screen_value = 42.392 - screen_variable * (7.796 - 1.961 * screen_variable);
  }

  return screen_value;
} 

Float ScreenFunction2(Float screen_variable) {

  // compute the value of the screening function 1.5*PHI1 - 0.5*PHI2

  Float screen_value;

  if (screen_variable > 1.0) {
    screen_value = 42.24 - 8.368 * log(screen_variable + 0.952);
  } else {
    screen_value = 41.734 - screen_variable * (6.484 - 1.250 * screen_variable);
  }

  return screen_value;
}

} // End namespace na63