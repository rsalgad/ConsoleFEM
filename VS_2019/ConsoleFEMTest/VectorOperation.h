#pragma once
#include "Matrix.h"
#include "Node.h"

class VectorOperation
{
public:
	static Matrix UnitVecX();
	static Matrix UnitVecY();
	static Matrix UnitVecZ();
	static Matrix RotateVectorX(Matrix &vec, double angle);
	static Matrix RotateVectorY(Matrix &vec, double angle);
	static Matrix RotateVectorZ(Matrix &vec, double angle);
	static Matrix RotateVectorIn3D(Matrix &vec, double angleX, double angleY, double angleZ);
	static Matrix VectorFromNodes(Node &n1, Node &n2);
	static Matrix VectorCrossProduct(Matrix &v1, Matrix &v2);
	static Matrix UnitVector(Matrix &vec);
	static double VectorLength(Matrix &vec);

private:
	VectorOperation();
};

