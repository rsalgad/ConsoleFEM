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
	static Matrix FullCholesky(Matrix &m, Matrix &f);
	static Matrix I(int dim);
	static Matrix Transpose(Matrix &m);
	static void GetCofactor(Matrix &matrix, Matrix &temp, int p, int q, int n);
	static Matrix GetInverseWithAdjoint(Matrix &m);
	static Matrix GetAdjointMatrix(Matrix &matrix);
	static double DeterminantWithCofactor(Matrix &m, int n);
	static double CalculateDeterminant(Matrix &m);
	static Matrix AddLineAndColumnWithTermAtPosition(Matrix &ori, int posLine, int posCol, double val);
	static double** CopyMatrixDouble(Matrix &m);
	static Matrix DeleteRowAndColumn(Matrix &m, int pos);
	static Matrix DeleteRow(Matrix &m, int pos);
	static Matrix CholeskyDecomposition(Matrix &m);
	static Matrix CholeskyDecompositionThreads(Matrix &m);
	static Matrix ForwardSubstitution(Matrix &L, Matrix &b);
	static Matrix BackSubstitution(Matrix &L, Matrix &b);
	static Matrix AddMatrixBottom(Matrix &ori, Matrix &toAdd);
	static Matrix AddMatrixRight(Matrix &ori, Matrix &toAdd);
	static void AddMatrixAtPosition(Matrix &ori, Matrix &toAdd, int initRowPos, int initColPos);
	static Matrix GetInverse(Matrix &m);
	static Matrix GetNormalizedUpperTriangular(Matrix &m);
	static Matrix GetNormalizedLowerTriangular(Matrix &m);
	static Matrix ExtractMatrixFromEnd(Matrix &m, int dim);
	static Matrix ExtractSetOfColumns(Matrix &ori, int nCols, int initCol, int step);
	static Matrix GetReducedMatrix(Matrix &GlobalMatrix, std::vector<Support> &vecSup, std::vector<Node> &vecNode);
	static double GetBiggestDiagTerm(Matrix& m);
	static void SwapLine(Matrix &matrix, int origIndex, int finalIndex);
	static void SwapColumn(Matrix &matrix, int origIndex, int finalIndex);
	static void MoveLineToEnd(Matrix &matrix, int origIndex);
	static void MoveColumnToEnd(Matrix &matrix, int origIndex);
	static void MoveLineAndColumnToEnd(Matrix &matrix, int origIndex);
	static Matrix ExtractMatrixBasedOnLineAndColumns(Matrix &matrix, int iniLine, int finalLine, int iniCol, int finalCol);
	static void ConvertToEigenMatrix(Matrix &m, Eigen::MatrixXd& eigenMatrix);
	static Matrix Sqrt(Matrix &m);
	static void PopulateDiagonalOnly(std::vector<double>& terms, Matrix &m);
	static Matrix FillMatrixBasedOnOtherMatrix(Matrix &toFill, Matrix &filler);

private:
	MatrixOperation();
};

