#include "pch.h"
#include "VectorOperation.h"

//<summary> Creates a unit vector pointing in the global x direction </summary> 
Matrix VectorOperation::UnitVecX() {
	Matrix vec = Matrix(3, 1);
	vec.GetMatrixDouble()[0][0] = 1;
	return vec;
}

//<summary> Creates a unit vector pointing in the global y direction </summary> 
Matrix VectorOperation::UnitVecY() {
	Matrix vec = Matrix(3, 1);
	vec.GetMatrixDouble()[1][0] = 1;
	return vec;
}

//<summary> Creates a unit vector pointing in the global z direction </summary> 
Matrix VectorOperation::UnitVecZ() {
	Matrix vec = Matrix(3, 1);
	vec.GetMatrixDouble()[2][0] = 1;
	return vec;
}

//<summary> Returns a vector that connects n2 and n1. The vector will point towards n2. </summary>
//<n1> The origin node of the vector </n1>
//<n2> The 'tip' node of the vector </n2>
Matrix VectorOperation::VectorFromNodes(const Node &n1, const Node &n2) {
	double** vector = Matrix::CreateMatrixDouble(3, 1);
	vector[0][0] = n2.GetX() - n1.GetX();
	vector[1][0] = n2.GetY() - n1.GetY();
	vector[2][0] = n2.GetZ() - n1.GetZ();
	return Matrix(vector, 3, 1);
}

//<summary> Performs the cross product of two vectors </summary>
Matrix VectorOperation::VectorCrossProduct(const Matrix &n1, const Matrix &n2) {
	double** vector = Matrix::CreateMatrixDouble(3, 1);
	vector[0][0] = n1.GetMatrixDouble()[1][0] * n2.GetMatrixDouble()[2][0] - n2.GetMatrixDouble()[1][0] * n1.GetMatrixDouble()[2][0];
	vector[1][0] = -(n1.GetMatrixDouble()[0][0] * n2.GetMatrixDouble()[2][0] - n2.GetMatrixDouble()[0][0] * n1.GetMatrixDouble()[2][0]);
	vector[2][0] = n1.GetMatrixDouble()[0][0] * n2.GetMatrixDouble()[1][0] - n2.GetMatrixDouble()[0][0] * n1.GetMatrixDouble()[1][0];
	return Matrix(vector, 3, 1);
}

//<summary> Calculates the length of a given vector </summary>
double VectorOperation::VectorLength(const  Matrix &vec) {
	return sqrt(pow(vec.GetMatrixDouble()[0][0], 2) + pow(vec.GetMatrixDouble()[1][0], 2) + pow(vec.GetMatrixDouble()[2][0], 2));
}

//<summary>Rotates a vector in the global X direction </summary>
//<angle>The angle to rotate, in radians </angle>
Matrix VectorOperation::RotateVectorX(const Matrix &vec, const double angle) {
	Matrix rot(3);
	double cos = std::cos(angle);
	double sin = std::sin(angle);
	rot.GetMatrixDouble()[0][0] = 1;
	rot.GetMatrixDouble()[1][1] = cos;
	rot.GetMatrixDouble()[1][2] = -sin;
	rot.GetMatrixDouble()[2][1] = sin;
	rot.GetMatrixDouble()[2][2] = cos;
	Matrix ans = rot * vec;

	return ans;
}

//<summary>Rotates a vector in the global Y direction </summary>
//<angle>The angle to rotate, in radians </angle>
Matrix VectorOperation::RotateVectorY(const Matrix &vec, const double angle) {
	Matrix rot(3);
	double cos = std::cos(angle);
	double sin = std::sin(angle);
	rot.GetMatrixDouble()[0][0] = cos;
	rot.GetMatrixDouble()[0][2] = sin;
	rot.GetMatrixDouble()[1][1] = 1;
	rot.GetMatrixDouble()[2][0] = -sin;
	rot.GetMatrixDouble()[2][2] = cos;
	Matrix ans = rot * vec;

	return ans;
}

//<summary>Rotates a vector in the global Z direction </summary>
//<angle>The angle to rotate, in radians </angle>
Matrix VectorOperation::RotateVectorZ(const Matrix &vec, const double angle) {
	Matrix rot(3);
	double cos = std::cos(angle);
	double sin = std::sin(angle);
	rot.GetMatrixDouble()[0][0] = cos;
	rot.GetMatrixDouble()[0][1] = -sin;
	rot.GetMatrixDouble()[1][0] = sin;
	rot.GetMatrixDouble()[1][1] = cos;
	rot.GetMatrixDouble()[2][2] = 1;
	Matrix ans = rot * vec;

	return ans;
}

//<summary>Rotates a vector in all the 3 global directions at once </summary>
//<angleX>The angle to rotate in the global X, in radians </angleX>
//<angleY>The angle to rotate in the global Y, in radians </angleY>
//<angleZ>The angle to rotate in the global Z, in radians </angleZ>
Matrix VectorOperation::RotateVectorIn3D(const Matrix &vec, const double angleX, const double angleY, const double angleZ) {
	Matrix rotX(3);
	Matrix rotY(3);
	Matrix rotZ(3);
	double cosx = std::cos(angleX);
	double sinx = std::sin(angleX);
	double cosy = std::cos(angleY);
	double siny = std::sin(angleY);
	double cosz = std::cos(angleZ);
	double sinz = std::sin(angleZ);
	rotX.GetMatrixDouble()[0][0] = 1;
	rotX.GetMatrixDouble()[1][1] = cosx;
	rotX.GetMatrixDouble()[1][2] = -sinx;
	rotX.GetMatrixDouble()[2][1] = sinx;
	rotX.GetMatrixDouble()[2][2] = cosx;
	rotY.GetMatrixDouble()[0][0] = cosy;
	rotY.GetMatrixDouble()[0][2] = siny;
	rotY.GetMatrixDouble()[1][1] = 1;
	rotY.GetMatrixDouble()[2][0] = -siny;
	rotY.GetMatrixDouble()[2][2] = cosy;
	rotZ.GetMatrixDouble()[0][0] = cosz;
	rotZ.GetMatrixDouble()[0][1] = -sinz;
	rotZ.GetMatrixDouble()[1][0] = sinz;
	rotZ.GetMatrixDouble()[1][1] = cosz;
	rotZ.GetMatrixDouble()[2][2] = 1;

	Matrix ans = rotZ * rotY * rotX * vec;

	return ans;
}

//<summary>Returns the unit vector of a specified vector </summary>
Matrix VectorOperation::UnitVector(const Matrix &vec) {
	double length = VectorLength(vec);
	double** vector = Matrix::CreateMatrixDouble(3, 1);
	vector[0][0] = vec.GetMatrixDouble()[0][0] / length;
	vector[1][0] = vec.GetMatrixDouble()[1][0] / length;
	vector[2][0] = vec.GetMatrixDouble()[2][0] / length;
	return Matrix(vector, 3, 1);
}
