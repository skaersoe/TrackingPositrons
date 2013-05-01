#include "Geometry/Box.hh"
#include <cmath>

namespace na63 {

Box::Box(const char* n, ThreeVector p, ThreeVector d) 
    : Volume(n,BOX) {

  // Initialize position and dimensions as well as the vectors describing the box.

  pos = p;
  dim = d;

  pos_vector[0] = pos[0];
  pos_vector[1] = pos[1];
  pos_vector[2] = pos[2];

  x_vector[0] = 0.5*dim[0];
  x_vector[1] = 0.0;
  x_vector[2] = 0.0;

  y_vector[0] = 0.0;
  y_vector[1] = 0.5*dim[1];
  y_vector[2] = 0.0;

  z_vector[0] = 0.0;
  z_vector[1] = 0.0;
  z_vector[2] = 0.5*dim[2];

  // Initialize constant elements of the rotation matrices.

  x_rotation.setZero();
  y_rotation.setZero();
  z_rotation.setZero();

  x_rotation(0,0) = 1;
  y_rotation(1,1) = 1;
  z_rotation(2,2) = 1;
}

bool Box::Inside(ThreeVector particle_position) const {

  Eigen::Vector3f particle_vector;

  particle_vector[0] = particle_position[0];
  particle_vector[1] = particle_position[1];
  particle_vector[2] = particle_position[2];

  Eigen::Vector3f x_top;
  Eigen::Vector3f x_bottom;
  Eigen::Vector3f x_negative;

  Eigen::Vector3f y_top;
  Eigen::Vector3f y_bottom;
  Eigen::Vector3f y_negative;

  Eigen::Vector3f z_top;
  Eigen::Vector3f z_bottom;
  Eigen::Vector3f z_negative;

  Eigen::Vector3f x_bottom_distance;
  Eigen::Vector3f x_top_distance;

  Eigen::Vector3f y_bottom_distance;
  Eigen::Vector3f y_top_distance;

  Eigen::Vector3f z_bottom_distance;
  Eigen::Vector3f z_top_distance;

  int x_top_sign;
  int x_bottom_sign;

  int y_top_sign;
  int y_bottom_sign;

  int z_top_sign;
  int z_bottom_sign;

  // Check if the particle is between the two edges normal to the x_vector.
  // if not, return false.

  x_top = pos_vector + x_vector;
  x_bottom = pos_vector + ( - x_vector);

  x_top_distance = particle_vector - x_top;

  if ((x_vector.dot(x_top_distance)) > 0)
    x_top_sign = 1;
  else
    x_top_sign = -1;

  x_bottom_distance = particle_vector - x_bottom;

  x_negative = -x_vector;

  if ((x_negative.dot(x_bottom_distance)) > 0)
    x_bottom_sign = 1;
  else
    x_bottom_sign = -1;

  if (x_top_sign == x_bottom_sign)
    return false;

  // Check if the particle is between the two edges normal to the y_vector.
  // if not, return false.

  y_top = pos_vector + y_vector;
  y_bottom = pos_vector + ( - y_vector);

  y_top_distance = particle_vector - y_top;

  if ((y_vector.dot(y_top_distance)) > 0)
    y_top_sign = 1;
  else
    y_top_sign = -1;

  y_bottom_distance = particle_vector - y_bottom;  

  y_negative = -y_vector;

  if ((y_negative.dot(y_bottom_distance)) > 0)
    y_bottom_sign = 1;
  else
    y_bottom_sign = -1;

  if (y_top_sign == y_bottom_sign)
    return false;



  // Check if the particle is between the two edges normal to the z_vector.
  // if not, return false.

  z_top = pos_vector + z_vector;
  z_bottom = pos_vector + ( - z_vector);

  z_top_distance = particle_vector - z_top;

  if ((z_vector.dot(z_top_distance)) > 0)
    z_top_sign = 1;
  else
    z_top_sign = -1;

  z_bottom_distance = particle_vector - z_bottom;  

  z_negative = -z_vector;

  if ((z_negative.dot(z_bottom_distance)) > 0)
    z_bottom_sign = 1;
  else
    z_bottom_sign = -1;

  if (z_top_sign == z_bottom_sign)
    return false;

  return true;

}

void Box::rotate(Float x_deg, Float y_deg, Float z_deg) {

  // Set the rotation matrices and execute the rotation on the box vectors.

  setRotation(x_deg, y_deg, z_deg);

  x_vector = total_rotation * x_vector;
  y_vector = total_rotation * y_vector;
  z_vector = total_rotation * z_vector;

  return;

}

void Box::setRotation(Float x_deg, Float y_deg, Float z_deg) {

  // Set the rotation matrices with the given rotations and calculate the total rotation matrix.

  x_rotation(1,1) = cos(x_deg);
  x_rotation(1,2) = -(sin(x_deg));
  x_rotation(2,1) = sin(x_deg);
  x_rotation(2,2) = cos(x_deg);

  y_rotation(0,0) = cos(y_deg);
  y_rotation(0,2) = sin(y_deg);
  y_rotation(2,0) = -(sin(y_deg));
  y_rotation(2,2) = cos(y_deg);

  z_rotation(0,0) = cos(z_deg);
  z_rotation(0,1) = -(sin(z_deg));
  z_rotation(1,0) = sin(z_deg);
  z_rotation(1,1) = cos(z_deg);

  total_rotation = x_rotation * y_rotation * z_rotation;

  return;
}

void Box::SetSpecificParameters(void *parameters) {
  BoxPars* box = (BoxPars*)parameters;
  pos.GPU(box->position);
  dim.GPU(box->dimensions);
  x_vector.GPU(box->x_vector);
  y_vector.GPU(box->y_vector);
  z_vector.GPU(box->z_vector);
}

InsideFunction inside_function() const {
  return NULL;
}

}
