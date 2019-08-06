#include "pch.h"
#include "Matrix.h"
#include <thread>
#include <mutex>
#include <vector>
#include <iostream>

Matrix::Matrix()
{
}

Matrix::~Matrix()
{
	if (!hasBeenCopied) {
		this->DestroyMatrixDouble();
	}
}

void Matrix::Print() {
	std::string matrix = this->ToString();
	std::cout << matrix << std::endl;
}

//Copy constructor
Matrix::Matrix(const Matrix& obj) {
	this->_matrix = obj.GetMatrixDouble();
	this->_dimX = obj._dimX;
	this->_dimY = obj._dimY;
	obj.hasBeenCopied = true;
}

//<summary> Initializes a matrix by passing its two dimensions </summary> 
//<dimX> Number of lines </dimX>
//<dimY> Number of columns </dimY>
Matrix::Matrix(int dimX, int dimY)
{
	_matrix = Matrix::CreateMatrixDouble(dimX, dimY);
	_dimX = dimX;
	_dimY = dimY;
}

//<summary> Initializes a square matrix by passing one dimension </summary> 
//<dim> Dimension </dim>
Matrix::Matrix(int dim)
{
	_matrix = Matrix::CreateMatrixDouble(dim);
	_dimX = dim;
	_dimY = dim;
}

//<summary> Initializes a matrix by passing its double array of contents and its two dimensions </summary> 
//<matrix> Previously created double** with the values of the matrix </matrix>
//<dimX> Number of lines </dimX>
//<dimY> Number of columns </dimY>
Matrix::Matrix(double** matrix, int dimX, int dimY) {
	_matrix = matrix;
	_dimX = dimX;
	_dimY = dimY;
}

//<summary> Initializes a square matrix by passing its double array of contents and its dimension </summary> 
//<matrix> Previously created double** with the values of the matrix </matrix>
//<dim> Dimension </dim>
Matrix::Matrix(double** matrix, int dim) {
	_matrix = matrix;
	_dimX = dim;
	_dimY = dim;
}

void Matrix::operator =(Matrix const &m2) {
	this->DestroyMatrixDouble();
	this->_matrix = m2.GetMatrixDouble();
	this->_dimX = m2.GetDimX();
	this->_dimY = m2.GetDimY();
	m2.hasBeenCopied = true;
}

void Matrix::SetDimensions(int dimX, int dimY) {
	_matrix = Matrix::CreateMatrixDouble(dimX, dimY);
	_dimX = dimX;
	_dimY = dimY;
}

//<summary> Handles how the matrix is transformed to strings when it needs to be printed somewhere </summary> 
std::string Matrix::ToString() {
	std::string row = "";
	if (_dimY != 1) {
		for (int i = 0; i < _dimX; i++) {
			for (int j = 0; j < _dimY; j++) {
				if (j == 0) {
					row += "[";
					row += std::to_string(_matrix[i][j]);
					row += " ";
				}
				else if (j == _dimY - 1) {
					row += std::to_string(_matrix[i][j]);
					row += "]\n";
				}
				else {
					row += std::to_string(_matrix[i][j]);
					row += " ";
				}
			}
		}
	}
	else {
		for (int i = 0; i < _dimX; i++) {
			row += "[";
			row += std::to_string(_matrix[i][0]);
			row += "]\n";
		}
	}
	return row;
}

//<summary> Returns the double** of the values of a matrix </summary> 
double** Matrix::GetMatrixDouble() const {
	return _matrix;
}

//<summary> Changes the double** of values of a matrix by the specified double** </summary> 
//<matrix> The double** that will replace the matrix original double** </matrix> 
void Matrix::SetMatrixDouble(double** matrix) {
	Matrix::DestroyMatrixDouble(_matrix, _dimX);
	_matrix = matrix;
}

//<summary> Returns the number of lines in a matrix </summary> 
int Matrix::GetDimX() const {
	return _dimX;
}

//<summary> Returns the number of columns in a matrix </summary> 
int Matrix::GetDimY() const {
	return _dimY;
}

//<summary> Creates a double** from the specified dimensions. </summary> 
//<dimX> Number of lines </dimX>
//<dimY> Number of columns </dimY>
double** Matrix::CreateMatrixDouble(int dimX, int dimY) {
	double** m = new double*[dimX];
	for (int i = 0; i < dimX; i++)
		m[i] = new double[dimY] {}; // '{}' initializes everything to zero
	return m;
}

//<summary> Creates a square double** from the specified dimension. </summary> 
//<dim> Dimension </dim>
double** Matrix::CreateMatrixDouble(int dim) {
	return CreateMatrixDouble(dim, dim);
}

//<summary> Destroys the double** of the matrix. </summary> 
//<comment> This is required because double** is created using the 'new' keyword, which stores it in the heap, which means it does not get deleted once it goes out of scope. If it is not deleted manually, there will be memory leaks. </comment> 
void Matrix::DestroyMatrixDouble() {
	for (int i = 0; i < _dimX; i++) {
		delete[] _matrix[i];
	}
	delete[] _matrix;
}

//<summary> Destroys the specified double**. </summary> 
//<dimX> The line dimension of the double** </dimX> 
void Matrix::DestroyMatrixDouble(double** d, int dimX) {
	for (int i = 0; i < dimX; i++) {
		delete[] d[i];
	}
	delete[] d;
}

//<summary> Calculate the sum of two matrices</summary>
Matrix Matrix::operator +(Matrix const &m2) {
	double** addedMatrix = Matrix::CreateMatrixDouble(_dimX, _dimY);
	for (int i = 0; i < _dimX; i++)
	{
		for (int j = 0; j < _dimY; j++)
		{
			addedMatrix[i][j] = _matrix[i][j] + m2._matrix[i][j];
		}
	}
	Matrix sum(addedMatrix, _dimX, _dimY);
	return sum;
}

//<summary> Calculate the '+=' of two matrices </summary>
void Matrix::operator +=(Matrix const &m2) {
	for (int i = 0; i < _dimX; i++)
	{
		for (int j = 0; j < _dimY; j++)
		{
			_matrix[i][j] += m2._matrix[i][j];
		}
	}
}

//<summary> Calculate the subtraction of two matrices </summary>
Matrix Matrix::operator -(Matrix const &m2) {
	double** subMatrix = Matrix::CreateMatrixDouble(_dimX, _dimY);
	for (int i = 0; i < _dimX; i++)
	{
		for (int j = 0; j < _dimY; j++)
		{
			subMatrix[i][j] = _matrix[i][j] - m2._matrix[i][j];
		}
	}
	Matrix sub(subMatrix, _dimX, _dimY);
	return sub;
}

//<summary> Calculate the multiplication of two matrices </summary>
Matrix Matrix::operator *(Matrix const &m2) {
	double** multMatrix = Matrix::CreateMatrixDouble(_dimX, m2._dimY);

	for (int i = 0; i < _dimX; i++)
	{
		for (int k = 0; k < m2._dimY; k++)
		{
			double multiplication = 0;
			for (int j = 0; j < _dimY; j++) {
				multiplication += _matrix[i][j] * m2._matrix[j][k];
			}
			multMatrix[i][k] = multiplication;
		}
	}
	Matrix mult(multMatrix, _dimX, m2._dimY);
	return mult;
}

//<summary> Calculate the multiplication of a matrix by a double </summary>
Matrix Matrix::operator *(double const &d) {
	double** multMatrix = Matrix::CreateMatrixDouble(_dimX, _dimY);
	for (int i = 0; i < _dimX; i++)
	{
		for (int j = 0; j < _dimY; j++)
		{
			multMatrix[i][j] = _matrix[i][j] * d;
		}
	}
	Matrix mult(multMatrix, _dimX, _dimY);
	return mult;
}

//<summary> Calculate the multiplication of a matrix by an integer </summary>
Matrix Matrix::operator *(int const &n) {
	double** multMatrix = Matrix::CreateMatrixDouble(_dimX, _dimY);
	for (int i = 0; i < _dimX; i++)
	{
		for (int j = 0; j < _dimY; j++)
		{
			multMatrix[i][j] = _matrix[i][j] * n;
		}
	}
	Matrix mult(multMatrix, _dimX, _dimY);
	return mult;
}


//Calculates the sparsity of the matrix
double Matrix::Sparsity() {
	double countNonZero = 0;
	double countZero = 0;
	for (int i = 0; i < _dimX; i++) {
		for (int j = 0; j < _dimY; j++) {
			if (_matrix[i][j] == 0) {
				countZero++;
			}
			else {
				countNonZero++;
			}
		}
	}

	return countZero / (_dimX * _dimY);
}