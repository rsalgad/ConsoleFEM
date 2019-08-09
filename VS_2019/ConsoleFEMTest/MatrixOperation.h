#pragma once
#include <Eigen/EigenValues>
#include "Node.h"
#include "Support.h"
#include "Matrix.h"

class MatrixOperation
{
public:
	static void Cholesky1(int ini, int k, double& sum, double** matrix);
	static void Cholesky2(int ini, int k, int i, double& sum, double** matrix);
	static Matrix FullCholesky(const Matrix &m, const Matrix &f);
	static Matrix I(int dim);
	static Matrix Transpose(const Matrix &m);
	static void GetCofactor(const Matrix &matrix, Matrix &temp, int p, int q, int n);
	static Matrix GetInverseWithAdjoint(const Matrix &m);
	static Matrix GetAdjointMatrix(const Matrix &matrix);
	static double DeterminantWithCofactor(const Matrix &m, int n);
	static double CalculateDeterminant(const Matrix &m);
	static Matrix AddLineAndColumnWithTermAtPosition(const Matrix &ori, const int posLine, const int posCol, const double val);
	static double** CopyMatrixDouble(const Matrix &m);
	static Matrix DeleteRowAndColumn(const Matrix &m, const int pos);
	static Matrix DeleteRow(const Matrix &m, const int pos);
	static Matrix CholeskyDecomposition(const Matrix &m);
	static Matrix CholeskyDecompositionThreads(const Matrix &m);
	static Matrix ForwardSubstitution(const Matrix &L, const Matrix &b);
	static Matrix BackSubstitution(const Matrix &L, const Matrix &b);
	static Matrix AddMatrixBottom(const Matrix &ori, const Matrix &toAdd);
	static Matrix AddMatrixRight(const Matrix &ori, const Matrix &toAdd);
	static void AddMatrixAtPosition(Matrix &ori, const Matrix &toAdd, const int initRowPos, const int initColPos);
	static Matrix GetInverse(const Matrix &m);
	static Matrix GetNormalizedUpperTriangular(const Matrix &m);
	static Matrix GetNormalizedLowerTriangular(const Matrix &m);
	static Matrix ExtractMatrixFromEnd(const Matrix &m, const int dim);
	static Matrix ExtractSetOfColumns(const Matrix &ori, const int nCols, const int initCol, const int step);
	static Matrix GetReducedMatrix(const Matrix &GlobalMatrix, std::vector<Support> &vecSup, std::vector<Node> &vecNode);
	static double GetBiggestDiagTerm(const Matrix& m);
	static void SwapLine(Matrix &matrix, const int origIndex, const int finalIndex);
	static void SwapColumn(Matrix &matrix, const int origIndex, const int finalIndex);
	static void MoveLineToEnd(Matrix &matrix, const int origIndex);
	static void MoveColumnToEnd(Matrix &matrix, const int origIndex);
	static void MoveLineAndColumnToEnd(Matrix &matrix, const int origIndex);
	static Matrix ExtractMatrixBasedOnLineAndColumns(const Matrix &matrix, const int iniLine, const int finalLine, const int iniCol, const int finalCol);
	static void ConvertToEigenMatrix(const Matrix &m, Eigen::MatrixXd& eigenMatrix);
	static Matrix Sqrt(const Matrix &m);
	static void PopulateDiagonalOnly(std::vector<double>& terms, Matrix &m);
	static Matrix FillMatrixBasedOnOtherMatrix(const Matrix &toFill, const Matrix &filler);

private:
	MatrixOperation();
};

