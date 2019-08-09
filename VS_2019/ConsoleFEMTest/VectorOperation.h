#pragma once
#include "Matrix.h"
#include "Node.h"

class VectorOperation
{
public:
	static Matrix UnitVecX();
	static Matrix UnitVecY();
	static Matrix UnitVecZ();
	static Matrix RotateVectorX(const Matrix &vec, const double angle);
	static Matrix RotateVectorY(const Matrix &vec, const double angle);
	static Matrix RotateVectorZ(const Matrix &vec, const double angle);
	static Matrix RotateVectorIn3D(const Matrix &vec, const double angleX, const double angleY, const double angleZ);
	static Matrix VectorFromNodes(const Node &n1, const Node &n2);
	static Matrix VectorCrossProduct(const Matrix &v1, const Matrix &v2);
	static Matrix UnitVector(const Matrix &vec);
	static double VectorLength(const Matrix &vec);

private:
	VectorOperation();
};

