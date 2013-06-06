#ifndef NA63_PARTICLE_OBJECT
#define NA63_PARTICLE_OBJECT

#include "Types.hh"

namespace na63 {

  class ParticleObject {

  private:

    point_3d position;
    point_3d velocity;
    int charge;
    double mass;
    double momentum;
    double tx;
    double ty;

  public:

    ParticleObject(double particleType, double momentumMean, double momentumSpread);

    ~ParticleObject();

    point_3d getPos() {return position;}
    point_3d getVel() {return velocity;}
    double getMomentum() {return momentum;}
    double getTx() {return tx;}
    double getTy() {return ty;}

    void setPos(double x, double y, double z);
    void setVel(double x, double y, double z);

  };

}

#endif