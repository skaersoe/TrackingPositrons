#ifndef NA63_GEOMETRY_MATERIAL_H
#define NA63_GEOMETRY_MATERIAL_H

namespace na63 {

  typedef struct {
    float density;
  } MaterialPars;

  class Material {

  public:
    Material(float density) {
      pars_.density = density;
    }
    ~Material();

    float density() { return pars_.density; }
    MaterialPars pars() { return pars_; }

  private:
    MaterialPars pars_;

  };

}

#endif /* NA63_GEOMETRY_MATERIAL_H */