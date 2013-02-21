#include <iostream>

#include "Particle.h"

int main(int argc,char *argv[]) {
  std::cout << "Starting..." << std::endl;
  Particle p(0.1,0.2,0.3);
  std::cout << "Created particle x=" << p.x << " y=" << p.y << " z=" << p.z << std::endl;
  std::cout << "...ending." << std::endl;
  return 0;
}