#ifndef NA63_GEOMETRY_BOX_H
#define NA63_GEOMETRY_BOX_H

#ifndef __CUDACC__
#include <eigen3/Eigen/Core>
#endif /* __CUDACC__ */
#include "Geometry/Volume.hh"

namespace na63 {

typedef struct {
  GPUThreeVector center;
  GPUThreeVector x_vector;
  GPUThreeVector y_vector;
  GPUThreeVector z_vector;
  // Pad to volume parameter size
  char padding[(int)(
    VOLUME_PARAMETER_SIZE
    - 4*sizeof(ThreeVector)
  )];
} BoxPars;

#ifdef RUNNING_CPP11
// Some extra safety to shield against human errors
static_assert(sizeof(BoxPars) == VOLUME_PARAMETER_SIZE,
    "Incorrect parameter size of class derived from Volume");
#endif

class Box : public Volume {

private:

  #ifndef __CUDACC__

	Eigen::Vector3f pos_vector; // The position of the box in vector form.

	Eigen::Vector3f x_vector; //
	Eigen::Vector3f y_vector; // Vectors pointing from center to edge in each dimension.
	Eigen::Vector3f z_vector; //

	Eigen::Matrix3f x_rotation; 	//
	Eigen::Matrix3f y_rotation; 	// Rotation matrices for each axis
	Eigen::Matrix3f z_rotation; 	// (see: http://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions)
	Eigen::Matrix3f total_rotation; // Total rotation (= x * y * z) see link above.

  #endif /* __CUDACC__ */

  static InsideFunction inside_function_;

protected:
  virtual void SetSpecificParameters(void *parameters);
  virtual InsideFunction inside_function() const { return inside_function_; }

public:
  Box(const char* n, ThreeVector position, ThreeVector dimensions);
  #ifdef RUNNING_CPP11
  Box(const char* n,
      Float x_pos, Float y_pos, Float z_pos,
      Float x_dim, Float y_dim, Float z_dim)
      : Box(n,ThreeVector(x_dim,y_dim,z_dim),ThreeVector(x_pos,y_pos,z_pos)) {}
  #endif

  #ifndef __CUDACC__

	Eigen::Vector3f getx() {return x_vector ;} //
	Eigen::Vector3f gety() {return y_vector ;} // Returns the box vectors.
	Eigen::Vector3f getz() {return z_vector ;} //

  #endif /* __CUDACC__ */

	virtual bool Inside(const FourVector& position) const; // Returns 1 if the given point is inside the box, 0 otherwise.

	void Rotate(Float x_deg, Float y_deg, Float z_deg); // Rotates the box around each axis with the given degrees for each axis.

	void SetRotation(Float x_deg, Float y_deg, Float z_deg); // Helper function for rotate(), sets the rotation matrices.

};

} // End namespace na63

#endif /* NA63_GEOMETRY_BOX_H */