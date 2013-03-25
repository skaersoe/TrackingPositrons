#include <Eigen/Core>

#ifndef NA63_BOX_H
#define NA63_BOX_H

namespace NA63 {

typedef struct {
	float x;
	float y;
	float z;
} point;

class Box {

private:
	point pos;
	point dim;

	Eigen::Vector3f x_vector;
	Eigen::Vector3f y_vector;
	Eigen::Vector3f z_vector;

	Eigen::Matrix3f x_rotation;
	Eigen::Matrix3f y_rotation;
	Eigen::Matrix3f z_rotation;
	Eigen::Matrix3f total_rotation;

public:
	Box(float x_dim, float y_dim, float z_dim,
			float x_pos, float y_pos, float z_pos);

	point getDimension() { return dim ;}

	point getPosition() {return pos ;}

	bool inside(point particle_position);

	void rotate(float x_deg, float y_deg, float z_deg);

	void setRotation(float x_deg, float y_deg, float z_deg);

};

}

#endif