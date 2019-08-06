#pragma once
#include <string>

class Matrix
{
public:
	Matrix();
	~Matrix();
	Matrix(int dimX, int dimY);
	Matrix(int dim);
	Matrix(double** matrix, int dimX, int dimY);
	Matrix(double** matrix, int dim);
	Matrix(const Matrix &obj);
	std::string ToString();
	Matrix operator +(Matrix const &m2);
	void operator +=(Matrix const &m2);
	Matrix operator -(Matrix const &m2);
	Matrix operator *(Matrix const &m2);
	Matrix operator *(double const &d);
	Matrix operator *(int const &n);
	void operator =(Matrix const &m2);
	double** GetMatrixDouble() const;
	void SetMatrixDouble(double** matrix);
	void DestroyMatrixDouble();
	int GetDimX() const;
	int GetDimY() const;
	void SetDimensions(int dimX, int dimY);
	void Print();
	double Sparsity();

	static double** CreateMatrixDouble(int dimX, int dimY);
	static double** CreateMatrixDouble(int dim);
	static void DestroyMatrixDouble(double** d, int dimX);

private:
	double** _matrix;
	int _dimX;
	int _dimY;
	mutable bool hasBeenCopied = false;
};

