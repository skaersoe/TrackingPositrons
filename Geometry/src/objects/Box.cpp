#include "Box.h"
#include <cmath>

namespace NA63 {

Box::Box(float x_dim, float y_dim, float z_dim,
    float x_pos, float y_pos, float z_pos) {

  pos.x = x_pos;
  pos.y = y_pos;
  pos.z = z_pos;

  dim.x = x_dim;
  dim.y = y_dim;
  dim.z = z_dim;

  x_vector[0] = 0.5*x_dim;
  x_vector[1] = 0.0;
  x_vector[2] = 0.0;

  y_vector[0] = 0.0;
  y_vector[1] = 0.5*y_dim;
  y_vector[2] = 0.0;

  z_vector[0] = 0.0;
  z_vector[1] = 0.0;
  z_vector[2] = 0.5*z_dim;

  x_rotation.setZero();
  y_rotation.setZero();
  z_rotation.setZero();

  x_rotation(0,0) = 1;
  y_rotation(1,1) = 1;
  z_rotation(2,2) = 1;
}

bool Box::inside(point particle_position) {

  bool x = particle_position.x > (pos.x - 0.5 * dim.x) &&
    particle_position.x < (pos.x + 0.5 * dim.x);

  bool y = particle_position.y > (pos.y - 0.5 * dim.y) &&
    particle_position.y < (pos.y + 0.5 * dim.y);

  bool z = particle_position.z > (pos.z - 0.5 * dim.z) &&
    particle_position.z < (pos.z + 0.5 * dim.z);

  bool inside = x && y && z;

  return inside;

}

void Box::rotate(float x_deg, float y_deg, float z_deg) {

  setRotation(x_deg, y_deg, z_deg);

  x_vector = total_rotation * x_vector;
  y_vector = total_rotation * y_vector;
  z_vector = total_rotation * z_vector;

  return;

}

void Box::setRotation(float x_deg, float y_deg, float z_deg) {

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

}
