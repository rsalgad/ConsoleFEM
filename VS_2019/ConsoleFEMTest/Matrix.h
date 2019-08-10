#pragma once
#include <string>

class Matrix
{
public:
	Matrix();
	~Matrix();
	Matrix(const int dimX, const int dimY);
	Matrix(const int dim);
	Matrix(double** matrix, const int dimX, const int dimY);
	Matrix(double** matrix, const int dim);
	Matrix(const const Matrix &obj);
	std::string ToString() const;
	Matrix operator +(Matrix const &m2) const;
	void operator +=(Matrix const &m2) const;
	Matrix operator -(Matrix const &m2) const;
	Matrix operator *(Matrix const &m2) const;
	Matrix operator *(double const &d) const;
	Matrix operator *(int const &n) const;
	void operator =(Matrix const &m2);
	double** GetMatrixDouble() const; //could make this const (read-only) and create another function like "ModifyMatrixDouble" that can actually change the values of double**
	void SetMatrixDouble(double** matrix);
	void DestroyMatrixDouble();
	int GetDimX() const;
	int GetDimY() const;
	void SetDimensions(const int dimX, const int dimY);
	void Print() const;
	double Sparsity() const;

	static double** CreateMatrixDouble(const int dimX, const int dimY);
	static double** CreateMatrixDouble(const int dim);
	static void DestroyMatrixDouble(double** d, const int dimX);

private:
	double** _matrix;
	int _dimX;
	int _dimY;
	mutable bool hasBeenCopied = false;
};

