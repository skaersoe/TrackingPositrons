#include "ParticleObject.hh"
#include "Types.hh"
#include "TRandom3.h"
#include <cmath>
#include <cassert>
#include <iostream>
#include <sys/time.h>

namespace na63 {

  const double pi = 3.14159265359;

  timeval time;

  int timeRetval = gettimeofday(&time, NULL);

  TRandom3 random = TRandom3(time.tv_usec);

  ParticleObject::ParticleObject(double particleType, double momentumMean, double momentumSpread) {

    // Only two particle types exist. Assert given type is either one.

    assert (particleType == ELECTRONIC || particleType == MUONIC);

    position.x = 0;
    position.y = 0;
    position.z = 0;

    // Total momentum is a gaussian distribution around 6 Gev.
    momentum = random.Gaus(momentumMean, momentumSpread);
    velocity.z = momentum;

    // Now the total momentum is distributed between x,y,z by a gaussian
    // Mean value is all momentum goes to z.

    double distributedMomentum = std::abs(random.Gaus(0.0,0.1) * momentum);

    velocity.z = std::sqrt(momentum * momentum - distributedMomentum * distributedMomentum);

    double randomAngle = random.Uniform(0.0,2.0*pi);

    velocity.x = distributedMomentum * std::cos(randomAngle);

    velocity.y = distributedMomentum * std::sin(randomAngle);

    // Set transverse momentum.

    tx = velocity.x / velocity.z;
    ty = velocity.y / velocity.z;

    // Set the charge. 50/50 for +1 or -1.

    if (random.Uniform(0.0,1.0) > 0.5) {
        charge = 1;
    }
    else {
        charge = -1;
    }

    // Set mass according to specified particle type.

    mass = particleType;

  }

  ParticleObject::~ParticleObject() {}

  void ParticleObject::setPos(double x, double y, double z) {
    position.x = x;
    position.y = y;
    position.z = z;
  }

  void ParticleObject::setVel(double x, double y, double z) {
    velocity.x = x;
    velocity.y = y;
    velocity.z = z;
  }

}