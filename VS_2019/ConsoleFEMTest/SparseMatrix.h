#pragma once
#include "Matrix.h"
#include <map>

class SparseMatrix
{
public:
	SparseMatrix();
	SparseMatrix(const int dimX, const int dimY);
	~SparseMatrix();
	std::map<int, std::map<int, double>>* GetMap();
	void SetSize(const int size);
	int SparseSize(const Matrix& m);
	int IndexOfFirstLine(const int line);
	int ElementsInLineXBeforeColY(const int lineIndex, const int col);
	double** CreateMatrixDouble(const int size);

	static SparseMatrix ConvertToSparseMatrix(const Matrix &m);
	static SparseMatrix ConvertToSymmetricSparseMatrix(const Matrix& m);

private:
	std::map<int, std::map<int, double>> _matrix;
	//double** _matrix;
	int _size;
	int _dimX;
	int _dimY;
};

