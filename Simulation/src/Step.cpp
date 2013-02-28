#include "Simulation/Step.h"

void apply_magnetic_field(float* p, float dt, float q) {
  const static float B[] = {0.1,0.1,0.1};

  float p_x = p[0];
  float p_y = p[1];
  float p_z = p[2];

  // Cross product
  p[0] = p_x + dt * (q * (p_y * B[2] - p_z * B[1]) );
  p[1] = p_y + dt * (q * (p_z * B[0] - p_x * B[2]) );
  p[2] = p_z + dt * (q * (p_x * B[1] - p_y * B[0]) );
}

void timestep(simple_particle_t* p, float dt) {
  timestep(p->r,p->p,p->m,p->q,dt);
}
void timestep(float* r, float* p, float m, float q, float dt) {
  float r_x = r[0];
  float r_y = r[1];
  float r_z = r[2];

  apply_magnetic_field(p,dt,q);

  r[0] = r[0] + dt * p[0] / m;
  r[1] = r[1] + dt * p[1] / m;
  r[2] = r[2] + dt * p[2] / m;
}