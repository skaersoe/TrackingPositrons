#include <Eigen/Core>

#ifndef NA63_BOX_H
#define NA63_BOX_H

namespace na63 {

// Define a point in 3-dimensions.

typedef struct {
	float x;
	float y;
	float z;
} point;

class Box {

private:
	point pos; // Position and dimension as 3-dimensional points.
	point dim; //

	float bsphere; // Radius of the bounding sphere.

	Eigen::Vector3f pos_vector; // The position of the box in vector form.

	Eigen::Vector3f x_vector; //
	Eigen::Vector3f y_vector; // Vectors pointing from center to edge in each dimension.
	Eigen::Vector3f z_vector; //

	Eigen::Matrix3f x_rotation; 	//
	Eigen::Matrix3f y_rotation; 	// Rotation matrices for each axis
	Eigen::Matrix3f z_rotation; 	// (see: http://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions)
	Eigen::Matrix3f total_rotation; // Total rotation (= x * y * z) see link above.

public:
	Box(float x_dim, float y_dim, float z_dim,
			float x_pos, float y_pos, float z_pos); // Construct a box with given size and position.

	~Box(); // Destructor.

	point getDimension() { return dim ;} // Return the size of the box, independent of coordinates (length, width, height).

	point getPosition() {return pos ;} // Return the position of the box center point.

	Eigen::Vector3f getx() {return x_vector ;} //
	Eigen::Vector3f gety() {return y_vector ;} // Returns the box vectors.
	Eigen::Vector3f getz() {return z_vector ;} //

	bool inside(point particle_position); // Returns 1 if the given point is inside the box, 0 otherwise.

	void rotate(float x_deg, float y_deg, float z_deg); // Rotates the box around each axis with the given degrees for each axis.

	void setRotation(float x_deg, float y_deg, float z_deg); // Helper function for rotate(), sets the rotation matrices.

};

}

#endif