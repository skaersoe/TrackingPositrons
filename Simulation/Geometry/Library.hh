#ifndef NA63_LIBRARY_H
#define NA63_LIBRARY_H

#include <iostream>
#include <cmath>
#ifndef __CUDACC__
#include <TRandom3.h>
#endif /* __CUDACC__ */

namespace na63 {

#ifndef __CUDACC__
#define RUNNING_CPP11
#endif /* __CUDACC__ */

// Define nullptr for <C++11
#ifndef nullptr
#define nullptr NULL
#endif

typedef float Float;
#ifndef Int_t
typedef int Int_t;
#endif

// All volumes must be defined here
typedef enum {
  SPHERE,
  BOX,
  SIMPLEBOX
} VolumeType;

template <class A, class B>
inline
void Copy3D(A& a, B& b) {
  b[0] = a[0];
  b[1] = a[1];
  b[2] = b[2];
}

inline Float ElementaryChargeToCoulomb(const Float e) {
  return 1.602176565e-19 * e;
}

inline Float JouleToGev(const Float J) {
  return 6.24150934e9 * J;
}

inline Float CartesianToSpherical_R(const Float x, const Float y, const Float z) {
  return sqrt(pow(x,2) + pow(y,2) + pow(z,2));
}

inline Float CartesianToSpherical_Theta(const Float x, const Float y, const Float z) {
  return acos(z/CartesianToSpherical_R(x,y,z));
}

inline Float CartesianToSpherical_Phi(const Float x, const Float y) {
  return atan(y/x);
}

inline Float Gamma(const Float& beta) {
  return 1/sqrt(1-pow(beta,2));
}

inline Float Gamma(const Float& energy, const Float& mass) {
  return energy / mass;
}

inline Float Beta(const Float& gamma) {
  return sqrt(1-1/pow(gamma,2));
}

typedef Float GPUThreeVector[3];
typedef Float GPUFourVector[4];

class FourVector;

class ThreeVector {

private:
  Float vector[3];

public:
  ThreeVector(const Float a, const Float b, const Float c) {
    vector[0] = a;
    vector[1] = b;
    vector[2] = c;
  }
  #ifdef RUNNING_CPP11
  ThreeVector() : ThreeVector(0,0,0) {}
  #endif

  Float& operator[] (const int i) {
    return vector[i];
  }

  Float operator[] (const int i) const {
    return vector[i];
  }

  ThreeVector& operator=(const FourVector& fv);

  friend std::ostream& operator<<(std::ostream& os, const ThreeVector& tv) {
    os << "(" << tv.vector[0]
       << "," << tv.vector[1]
       << "," << tv.vector[2] << ")";
    return os;
  }

  friend ThreeVector operator-(const ThreeVector& lhs, const ThreeVector& rhs) {
    return ThreeVector(
             lhs.vector[0] - rhs.vector[0],
             lhs.vector[1] - rhs.vector[1],
             lhs.vector[2] - rhs.vector[2]
           );
  }

  ThreeVector& operator-=(const ThreeVector& rhs) {
    vector[0] -= rhs.vector[0];
    vector[1] -= rhs.vector[1];
    vector[2] -= rhs.vector[2];
    return *this;
  }

  Float length() const {
    return sqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2]);
  }

  void Normalize() {
    Float l = length();
    vector[0] /= l;
    vector[1] /= l;
    vector[2] /= l;
  }

  ThreeVector Normalized() const {
    Float l = length();
    return ThreeVector(vector[0]/l,vector[1]/l,vector[2]/l);
  }

  void Extend(Float a) {
    vector[0] *= a;
    vector[1] *= a;
    vector[2] *= a;
  }

  void Set(Float a, Float b, Float c) {
    vector[0] = a;
    vector[1] = b;
    vector[2] = c;
  }

  void Rotate(const ThreeVector& unit_vector) {

    Float u1 = unit_vector[0];
    Float u2 = unit_vector[1];
    Float u3 = unit_vector[2];
    Float up = u1*u1 + u2*u2;

    if (up) {
      up = sqrt(up);
      Float px = vector[0], py = vector[1], pz = vector[2];
      vector[0] = (u1*u3*px - u2*py + u1*up*pz)/up;
      vector[1] = (u2*u3*px + u1*py + u2*up*pz)/up;
      vector[2] = (u3*u3*px -    px + u3*up*pz)/up;
    } else {
      if (u3 < 0.0) {
        vector[0] = -vector[0]; vector[2] = -vector[2];
      }
    }

}

  // Insert into GPU format (array)
  void GPU(GPUThreeVector tv) const {
    tv[0] = vector[0];
    tv[1] = vector[1];
    tv[2] = vector[2];
  }

};

class FourVector {

private:
  Float vector[4];

public:

  Float& operator[] (const int i) {
    return vector[i];
  }

  Float operator[] (const int i) const {
    return vector[i];
  }

  FourVector(const Float a, const Float b, const Float c, const Float d) {
    vector[0] = a;
    vector[1] = b;
    vector[2] = c;
    vector[3] = d;
  }
  #ifdef RUNNING_CPP11
  FourVector() : FourVector(0,0,0,0) {}
  FourVector(ThreeVector tv) : FourVector(tv,0) {}
  FourVector(ThreeVector tv, Float f) : FourVector(tv[0],tv[1],tv[2],f) {}
  #endif

  Float R() const {
    return CartesianToSpherical_R(vector[0],vector[1],vector[2]);
  }
  Float Theta() const {
    return CartesianToSpherical_Theta(vector[0],vector[1],vector[2]);
  }
  Float Phi() const {
    return CartesianToSpherical_Phi(vector[0],vector[1]);
  }

  FourVector& operator=(const GPUFourVector& gpu_fv) {
    vector[0] = gpu_fv[0];
    vector[1] = gpu_fv[1];
    vector[2] = gpu_fv[2];
    vector[3] = gpu_fv[3];
    return *this;
  }

  FourVector& operator=(const ThreeVector& tv) {
    vector[0] = tv[0];
    vector[1] = tv[1];
    vector[2] = tv[2];
    vector[3] = 0;
    return *this;
  }

  friend std::ostream& operator<<(std::ostream& os, const FourVector& fv) {
    os << "(" << fv.vector[0]
       << "," << fv.vector[1]
       << "," << fv.vector[2]
       << "," << fv.vector[3] << ")";
    return os;
  }

  friend FourVector operator-(const FourVector& lhs, const FourVector& rhs) {
    return FourVector(
             lhs.vector[0] - rhs.vector[0],
             lhs.vector[1] - rhs.vector[1],
             lhs.vector[2] - rhs.vector[2],
             lhs.vector[3] - rhs.vector[3]
           );
  }

  FourVector& operator-=(const FourVector& rhs) {
    vector[0] -= rhs.vector[0];
    vector[1] -= rhs.vector[1];
    vector[2] -= rhs.vector[2];
    vector[3] -= rhs.vector[3];
    return *this;
  }

  inline void Set(Float a, Float b, Float c, Float d) {
    vector[0] = a;
    vector[1] = b;
    vector[2] = c;
    vector[3] = d;
  }

  inline void Normalize() {
    Float l = length();
    vector[0] /= l;
    vector[1] /= l;
    vector[2] /= l;
  }

  ThreeVector Normalized() const {
    Float l = length();
    return ThreeVector(vector[0]/l,vector[1]/l,vector[2]/l);
  }

  Float length() const;

  inline void Extend(Float a) {
    vector[0] *= a;
    vector[1] *= a;
    vector[2] *= a;
  }

  // Insert into GPU format (array)
  void GPU(GPUFourVector fv) const {
    fv[0] = vector[0];
    fv[1] = vector[1];
    fv[2] = vector[2];
    fv[3] = vector[3];
  }
};

inline ThreeVector SphericalToCartesian(
    const Float r, const Float theta, const Float phi) {
  return ThreeVector(r * sin(theta) * cos(phi), // x
                     r * sin(theta) * sin(phi), // y
                     r * cos(theta));           // z
}

class Geometry;
class Volume;
class Material;

typedef bool (*InsideFunction)(const GPUFourVector&,const void*);

template <class VectorType, class ArrayType>
void ParameterVectorToArray(VectorType *vec, ArrayType **arr) {
  int size = vec->size();
  *arr = new ArrayType[size];
  for (int i=0;i<size;i++) {
    (*arr)[i] = (*vec)[i].GPU();
  }
}

#define VOLUME_PARAMETER_SIZE 128 - 3*sizeof(int)
typedef struct {
  // Generic fields
  int function_index;
  int material_index;
  // Volume-specific fields. Derived classes must pad parameters to this size
  char specific[VOLUME_PARAMETER_SIZE];
} VolumePars;

} // End namespace na63

#endif /* NA63_LIBRARY_H */