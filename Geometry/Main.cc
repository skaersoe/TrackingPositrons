#include <iostream>
#include "Geometry/Geometry.hh"
#include "Geometry/Sphere.hh"

using namespace na63;

int main(void) {
  Geometry *geometry = new Geometry();
  Material *m_a = new Material(0.5);
  Material *m_b = new Material(0.01);
  geometry->AddMaterial(*m_a);
  geometry->AddMaterial(*m_b);
  Volume *v_a = new Sphere(m_a,{0,0,1},0.5);
  geometry->AddVolume(v_a);
  return 0;
}