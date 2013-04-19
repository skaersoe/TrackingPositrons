#ifndef NA63_GEOMETRY_MATERIAL_H
#define NA63_GEOMETRY_MATERIAL_H

#include <string>

namespace na63 {

  typedef struct {
    float density;
  } MaterialPars;

  class Material {

  public:
    Material(const char* n, float density) : name_(n) {
      pars_.density = density;
    }
    ~Material() {}

    float density() const { return pars_.density; }
    std::string name() { return name_; }
    MaterialPars pars() const { return pars_; }

  private:
    MaterialPars pars_;
    std::string name_;

  };

}

#endif /* NA63_GEOMETRY_MATERIAL_H */