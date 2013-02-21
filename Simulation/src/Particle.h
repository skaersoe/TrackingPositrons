class Particle {
  public:
  float x, y, z;
  float px, py, pz;
  float c, m;
  Particle(float, float, float);
};

Particle::Particle(float x, float y, float z) {
  this->x = x; this->y = y; this->z = z;
}