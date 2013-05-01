#ifndef NA63_GEOMETRY_BOX_H
#define NA63_GEOMETRY_BOX_H

#include <eigen3/Eigen/Core>
#include "Geometry/Volume.hh"

namespace na63 {

typedef struct {
	GPUThreeVector position;
	GPUThreeVector dimensions;
	GPUThreeVector x_vector;
	GPUThreeVector y_vector;
	GPUThreeVector z_vector;
  // Pad to volume parameter size
  char padding[(int)(
    VOLUME_PARAMETER_SIZE
    - 5*sizeof(GPUThreeVector)
  )];
} BoxPars;


#ifdef RUNNING_CPP11
// Some extra safety to shield against human errors
static_assert(sizeof(BoxPars) == VOLUME_PARAMETER_SIZE,
    "Incorrect parameter size of class derived from Volume");
#endif

// Define a ThreeVector in 3-dimensions.

class Box : public Volume {

private:
	ThreeVector pos; // Position and dimension as 3-dimensional ThreeVectors.
	ThreeVector dim; //

	// Float bsphere; // Radius of the bounding sphere.

	Eigen::Vector3f pos_vector; // The position of the box in vector form.

	Eigen::Vector3f x_vector; //
	Eigen::Vector3f y_vector; // Vectors ThreeVectoring from center to edge in each dimension.
	Eigen::Vector3f z_vector; //

	Eigen::Matrix3f x_rotation; 	//
	Eigen::Matrix3f y_rotation; 	// Rotation matrices for each axis
	Eigen::Matrix3f z_rotation; 	// (see: http://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions)
	Eigen::Matrix3f total_rotation; // Total rotation (= x * y * z) see link above.

public:
	Box(ThreeVector pos, ThreeVector dim);
	#ifdef RUNNING_CPP_11
	Box(Float x_dim, Float y_dim, Float z_dim,
			Float x_pos, Float y_pos, Float z_pos)
			: Box(ThreeVector(x_pos,y_pos,z_pos),ThreeVector(x_dim,y_dim,z_dim)) {}
	#endif

	~Box(); // Destructor.

	ThreeVector getDimension() { return dim ;} // Return the size of the box, independent of coordinates (length, width, height).

	ThreeVector getPosition() {return pos ;} // Return the position of the box center ThreeVector.

	Eigen::Vector3f getx() {return x_vector ;} //
	Eigen::Vector3f gety() {return y_vector ;} // Returns the box vectors.
	Eigen::Vector3f getz() {return z_vector ;} //

	virtual bool Inside(ThreeVector particle_position) const; // Returns 1 if the given ThreeVector is inside the box, 0 otherwise.

	virtual InsideFunction inside_function() const;

	void rotate(Float x_deg, Float y_deg, Float z_deg); // Rotates the box around each axis with the given degrees for each axis.

	void setRotation(Float x_deg, Float y_deg, Float z_deg); // Helper function for rotate(), sets the rotation matrices.

};

}

#endif /* NA63_GEOMETRY_BOX_H */