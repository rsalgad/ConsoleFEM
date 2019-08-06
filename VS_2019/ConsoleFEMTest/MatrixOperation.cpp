#include "pch.h"
#include "Matrix.h"
#include "MatrixOperation.h"
#include <iostream>
#include <mutex>
#include <thread>
//<summary> Creates an identity matrix with the desired dimension </summary> 
//<dim> Dimension </dim>
Matrix MatrixOperation::I(int dim) {

	double** m = Matrix::CreateMatrixDouble(dim);

	for (int i = 0; i < dim; i++) {
		m[i][i] = 1;
	}
	Matrix myMatrix(m, dim);
	return myMatrix;
}


//<summary> Copies the double** of a matrix and return it. </summary> 
//<comment> This is required because double** is a pointer, and changing its values before copying it will change all the matrices that uses that double** </comment> 
double** MatrixOperation::CopyMatrixDouble(Matrix &m) {
	double** matrix = Matrix::CreateMatrixDouble(m.GetDimX(), m.GetDimY());
	for (int i = 0; i < m.GetDimX(); i++) {
		for (int j = 0; j < m.GetDimY(); j++) {
			matrix[i][j] = m.GetMatrixDouble()[i][j];
		}
	}
	return matrix;
}

//<summary> Calculate the cofactor of a matrix. </summary>
//<comment> This function is part of a group of functions to calculate determinant and inverse of matrices. </comment>
void MatrixOperation::GetCofactor(Matrix &matrix, Matrix &temp, int p, int q, int n) {
	int i = 0, j = 0;

	// Looping for each element of the matrix 
	for (int row = 0; row < n; row++)
	{
		for (int col = 0; col < n; col++)
		{
			//  Copying into temporary matrix only those element 
			//  which are not in given row and column 
			if (row != p && col != q)
			{
				temp.GetMatrixDouble()[i][j++] = matrix.GetMatrixDouble()[row][col];

				// Row is filled, so increase row index and 
				// reset col index 
				if (j == n - 1)
				{
					j = 0;
					i++;
				}
			}
		}
	}
}

//<summary> Calculate the determinant of a matrix using the cofactor method. </summary>
//<comment> This function is part of a group of functions to calculate determinant and inverse of matrices. </comment>
double MatrixOperation::DeterminantWithCofactor(Matrix &m, int n) {
	double D = 0; // Initialize result 

	//  Base case : if matrix contains single element 
	if (n == 1)
		return m.GetMatrixDouble()[0][0];

	Matrix temp(m.GetDimX(), m.GetDimY()); // To store cofactors 

	int sign = 1;  // To store sign multiplier 

	// Iterate for each element of first row 
	for (int f = 0; f < n; f++)
	{
		// Getting Cofactor of m[0][f] 
		GetCofactor(m, temp, 0, f, n);
		D += sign * m.GetMatrixDouble()[0][f] * DeterminantWithCofactor(temp, n - 1);

		// terms are to be added with alternate sign 
		sign = -sign;
	}

	return D;
}

//<summary> Calculate the adjoint matrix of the specified matrix. </summary>
//<matrix> The matrix to calculate the adjoint matrix </matrix>
//<comment> This function is part of a group of functions to calculate determinant and inverse of matrices. </comment>
Matrix MatrixOperation::GetAdjointMatrix(Matrix &matrix) {
	int dim = matrix.GetDimX();
	Matrix adj(dim, dim);

	if (dim == 1)
	{
		adj.GetMatrixDouble()[0][0] = 1;
		return adj;
	}

	// temp is used to store cofactors of A[][] 
	int sign = 1;
	Matrix temp(dim, dim);

	for (int i = 0; i < dim; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			// Get cofactor of A[i][j] 
			GetCofactor(matrix, temp, i, j, dim);

			// sign of adj[j][i] positive if sum of row 
			// and column indexes is even. 
			sign = ((i + j) % 2 == 0) ? 1 : -1;

			// Interchanging rows and columns to get the 
			// transpose of the cofactor matrix 
			adj.GetMatrixDouble()[j][i] = (sign)*(DeterminantWithCofactor(temp, dim - 1));
		}
	}

	return adj;
}

//<summary> Calculate the determinant of a matrix </summary>
//<m> The matrix to calculate the determinant </m>
//<comment> This function is deprecated, use the CalculateDeterminantWithCofactor. This function can result in errors. </comment>
double MatrixOperation::CalculateDeterminant(Matrix &m) {
	double** matrix = Matrix::CreateMatrixDouble(m.GetDimX(), m.GetDimY());
	matrix = CopyMatrixDouble(m);
	for (int k = 0; k < m.GetDimX(); k++) // This index keeps zeroeing everything related to this row.
	{
		if (k < m.GetDimX() - 1) {
			for (int i = 0; i < m.GetDimX(); i++) {
				if (i > k) {
					double n1 = matrix[i][k] / matrix[k][k];
					for (int j = 0; j < m.GetDimY(); j++) {
						matrix[i][j] = matrix[i][j] - n1 * matrix[k][j];
					}
				}
			}
		}
	}
	double determinant = matrix[0][0];
	for (int i = 1; i < m.GetDimX(); i++) {
		determinant *= matrix[i][i];
	}
	return determinant;
}

//<summary> Calculate the inverse of a matrix using the adjoint matrix method </summary>
//<m> The matrix to calculate the inverse </m>
Matrix MatrixOperation::GetInverseWithAdjoint(Matrix &m) {
	int dim = m.GetDimX();
	Matrix inverse(dim, dim);

	// Find determinant of m[][] 
	double det = DeterminantWithCofactor(m, dim);
	if (det == 0)
	{
		std::cout << "Singular matrix, can't find its inverse" << std::endl;
		return false;
	}

	// Find adjoint 
	Matrix adj = GetAdjointMatrix(m);

	// Find Inverse using formula "inverse(A) = adj(A)/det(A)" 
	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			inverse.GetMatrixDouble()[i][j] = adj.GetMatrixDouble()[i][j] / det;
		}
	}

	return inverse;
}

//<summary> Calculate the transpose of a a matrix </summary>
//<m> The matrix to be transposed </m>
Matrix MatrixOperation::Transpose(Matrix &m) {
	double** transp = Matrix::CreateMatrixDouble(m.GetDimY(), m.GetDimX());
	for (int i = 0; i < m.GetDimX(); i++) {
		for (int j = 0; j < m.GetDimY(); j++) {
			transp[j][i] = m.GetMatrixDouble()[i][j];
		}
	}
	Matrix transpose(transp, m.GetDimY(), m.GetDimX());
	return transpose;
}

//<summary> Calculate theinverse of a a matrix</summary>
//<m> The matrix to calculate the inverse </m>
//<comment> This function is deprecated, use the GetInverseWithAdjoint. This function can result in errors. </comment>
Matrix MatrixOperation::GetInverse(Matrix &m) {
	Matrix In = I(m.GetDimX());
	Matrix expanded = AddMatrixRight(m, In);
	Matrix lower = GetNormalizedLowerTriangular(expanded);
	Matrix upper = GetNormalizedUpperTriangular(lower);
	Matrix finalMatrix = ExtractMatrixFromEnd(upper, m.GetDimX());

	return finalMatrix;
}



//<summary> Deletes a row and a column of the specified matrix </summary>
//<m> The specified matrix </m>
//<pos> The position (in array notation) to delete</pos>
Matrix MatrixOperation::DeleteRowAndColumn(Matrix &m, int pos) {
	int sizeX = m.GetDimX() - 1;
	int sizeY = m.GetDimY() - 1;
	double** reduced = Matrix::CreateMatrixDouble(sizeX, sizeY);
	for (int i = 0; i < m.GetDimX(); i++) {
		for (int j = 0; j < m.GetDimY(); j++) {
			if (i < pos && j < pos) {
				reduced[i][j] = m.GetMatrixDouble()[i][j];
			}
			else if (i == pos || j == pos) {
				//Do Nothing
			}
			else if (i < pos && j >= pos) {
				reduced[i][j - 1] = m.GetMatrixDouble()[i][j];
			}
			else {
				if (j < pos) {
					reduced[i - 1][j] = m.GetMatrixDouble()[i][j];
				}
				else {
					reduced[i - 1][j - 1] = m.GetMatrixDouble()[i][j];
				}
			}
		}
	}
	return Matrix(reduced, sizeX, sizeY);
}

//<summary> Deletes a row of the specified matrix </summary>
//<m> The specified matrix </m>
//<pos> The position (in array notation) to delete</pos>
Matrix MatrixOperation::DeleteRow(Matrix &m, int pos) {
	int sizeX = m.GetDimX() - 1; //will delete 1 row
	int sizeY = m.GetDimY(); //won't change columns
	double** reduced = Matrix::CreateMatrixDouble(sizeX, sizeY);
	for (int i = 0; i < m.GetDimX(); i++) {
		for (int j = 0; j < m.GetDimY(); j++) {
			if (i < pos) {
				reduced[i][j] = m.GetMatrixDouble()[i][j];
			}
			else if (i == pos) {
				//Do Nothing
			}
			else {
				reduced[i - 1][j] = m.GetMatrixDouble()[i][j];
			}
		}
	}
	return Matrix(reduced, sizeX, sizeY);
}

//<summary> Performs the Cholesky decomposition of a matrix and returns the L matrix, where LL'=A </summary>
//<m> The matrix to calculate the cholesky</m>
Matrix MatrixOperation::CholeskyDecomposition(Matrix &m) {
	int sizeX = m.GetDimX();
	int sizeY = m.GetDimY();
	double** matrix = Matrix::CreateMatrixDouble(sizeX, sizeY);

	for (int i = 0; i < sizeX; i++)
	{
		for (int k = 0; k <= i; k++)
		{
			if (k == i)
			{
				double sum = 0;
				for (int j = 0; j < k; j++)
				{
					sum += matrix[k][j] * matrix[k][j];
				}
				matrix[k][k] = sqrt(m.GetMatrixDouble()[k][k] - sum);
			}
			else
			{
				double sum = 0;
				for (int j = 0; j < k; j++)
				{
					sum += matrix[i][j] * matrix[k][j];
				}
				matrix[i][k] = 1.0 / matrix[k][k] * (m.GetMatrixDouble()[i][k] - sum);
			}
		}
	}
	return Matrix(matrix, sizeX, sizeY);
}

//<summary> Performs the Cholesky decomposition of a matrix and returns the L matrix, where LL'=A </summary>
//<m> The matrix to calculate the cholesky</m>
//<comment> This function does not work as it is now. Still needs better implementation</comment>
Matrix MatrixOperation::CholeskyDecompositionThreads(Matrix &m) {
	int sizeX = m.GetDimX();
	int sizeY = m.GetDimY();
	double** matrix1 = Matrix::CreateMatrixDouble(sizeX, sizeY);
	double** matrix2 = Matrix::CreateMatrixDouble(sizeX, sizeY);
	double** matrix3 = Matrix::CreateMatrixDouble(sizeX, sizeY);
	double** matrix4 = Matrix::CreateMatrixDouble(sizeX, sizeY);

	for (int i = 0; i < sizeX; i++)
	{
		for (int k = 0; k <= i; k++)
		{
			if (k == i)
			{
				double sum = 0;
				if (k > 3) {
					double nThreads = 4.0;
					double amount = k / nThreads;

					int size1, size2;
					if (fmod(k, nThreads) != 0) {
						size1 = floor(amount);
						size2 = fmod(amount, nThreads) * nThreads;
					}
					else {
						size1 = amount;
						size2 = 0;
					}


					std::mutex matrixLock;
					std::vector<std::thread> threadList;
					double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0;
					threadList.emplace_back(Cholesky1, 0, size1, std::ref(sum1), std::ref(matrix1));
					threadList.emplace_back(Cholesky1, size1, 2 * size1, std::ref(sum1), std::ref(matrix2));
					threadList.emplace_back(Cholesky1, 2 * size1, 3 * size1, std::ref(sum1), std::ref(matrix3));
					Cholesky1(3 * size1, k, sum4, matrix4);

					for (int i = 0; i < threadList.size(); i++) {
						threadList[i].join();
					}
					sum = sum1 + sum2 + sum3 + sum4;
				}
				else {
					for (int j = 0; j < k; j++)
					{
						sum += matrix1[k][j] * matrix1[k][j];
					}
				}
				matrix1[k][k] = sqrt(m.GetMatrixDouble()[k][k] - sum);
			}
			else
			{
				double sum = 0;
				if (k > 3) {
					double nThreads = 4.0;
					double amount = k / nThreads;

					int size1, size2;
					if (fmod(k, nThreads) != 0) {
						size1 = floor(amount);
						size2 = fmod(amount, nThreads) * nThreads;
					}
					else {
						size1 = amount;
						size2 = 0;
					}


					std::mutex matrixLock;
					std::vector<std::thread> threadList;
					double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0;
					threadList.emplace_back(Cholesky2, 0, size1, i, std::ref(sum1), std::ref(matrix1));
					threadList.emplace_back(Cholesky2, size1, 2 * size1, i, std::ref(sum1), std::ref(matrix2));
					threadList.emplace_back(Cholesky2, 2 * size1, 3 * size1, i, std::ref(sum1), std::ref(matrix3));
					Cholesky2(3 * size1, k, i, sum4, matrix4);

					for (int i = 0; i < threadList.size(); i++) {
						threadList[i].join();
					}
					sum = sum1 + sum2 + sum3 + sum4;
				}
				else {
					for (int j = 0; j < k; j++)
					{
						sum += matrix1[i][j] * matrix1[k][j];
					}
				}
				matrix1[i][k] = 1.0 / matrix1[k][k] * (m.GetMatrixDouble()[i][k] - sum);
			}
		}
	}
	return Matrix(matrix1, sizeX, sizeY);
}

//<summary> Helps calculate the Cholesky decomposition with threads, see CholeskyDecompositionThreads </summary>
void MatrixOperation::Cholesky1(int ini, int k, double& sum, double** matrix) {
	for (int j = ini; j < k; j++)
	{
		sum += matrix[k][j] * matrix[k][j];
	}
}

//<summary> Helps calculate the Cholesky decomposition with threads, see CholeskyDecompositionThreads </summary>
void MatrixOperation::Cholesky2(int ini, int k, int i, double& sum, double** matrix) {
	for (int j = ini; j < k; j++)
	{
		sum += matrix[i][j] * matrix[k][j];
	}
}

//<summary> Performs forward substitution matrix operations to solve a system of the type Lx=b. Returns x. </summary>
Matrix MatrixOperation::ForwardSubstitution(Matrix &L, Matrix &b) {
	//performs the Forward substitution on a system of the type Lx=b, and returns the x matrix
	//L needs to be a lower triangle
	int size = b.GetDimX();
	double** matrix = Matrix::CreateMatrixDouble(size, 1);

	for (int i = 0; i < size; i++)
	{
		double sum = 0;
		for (int j = 0; j < i; j++)
		{
			sum += L.GetMatrixDouble()[i][j] * matrix[j][0];
		}
		matrix[i][0] = (b.GetMatrixDouble()[i][0] - sum) / L.GetMatrixDouble()[i][i];
	}
	return Matrix(matrix, size, 1);
}

//<summary> Performs back substitution matrix operations to solve a system of the type Lx=b. Returns x. </summary>
Matrix MatrixOperation::BackSubstitution(Matrix &L, Matrix &b) {
	//performs the Back substitution on a system of the type Lx=b, and returns the x matrix
	//L needs to be a lower triangle
	int size = b.GetDimX();
	double** matrix = Matrix::CreateMatrixDouble(size, 1);

	for (int i = size - 1; i >= 0; i--)
	{
		double sum = 0;
		for (int j = size - 1; j > i; j--)
		{
			sum += L.GetMatrixDouble()[i][j] * matrix[j][0];
		}
		matrix[i][0] = (b.GetMatrixDouble()[i][0] - sum) / L.GetMatrixDouble()[i][i];
	}
	return Matrix(matrix, size, 1);
}

//<summary> Performs all the steps in the Cholesky decomposition on the system F=Kd. Returns d. </summary>
//<m> The K matrix </m>
//<f> The F matrix </f>
Matrix MatrixOperation::FullCholesky(Matrix &m, Matrix &f) {
	int sizeX = m.GetDimX();
	int sizeY = m.GetDimY();
	double** matrix = Matrix::CreateMatrixDouble(sizeX, sizeY);

	for (int i = 0; i < sizeX; i++)
	{
		for (int k = 0; k <= i; k++)
		{
			if (k == i)
			{
				double sum = 0;
				for (int j = 0; j < k; j++)
				{
					sum += matrix[k][j] * matrix[k][j];
				}
				matrix[k][k] = sqrt(m.GetMatrixDouble()[k][k] - sum);
			}
			else
			{
				double sum = 0;
				for (int j = 0; j < k; j++)
				{
					sum += matrix[i][j] * matrix[k][j];
				}
				matrix[i][k] = 1.0 / matrix[k][k] * (m.GetMatrixDouble()[i][k] - sum);
			}
		}
	}

	Matrix cholesky(matrix, sizeX, sizeY);
	Matrix choleskyTransp = Transpose(cholesky);

	int size = f.GetDimX();
	double** forwardSub = Matrix::CreateMatrixDouble(f.GetDimX(), 1);

	for (int i = 0; i < size; i++)
	{
		double sum = 0;
		for (int j = 0; j < i; j++)
		{
			sum += cholesky.GetMatrixDouble()[i][j] * forwardSub[j][0];
		}
		forwardSub[i][0] = (f.GetMatrixDouble()[i][0] - sum) / cholesky.GetMatrixDouble()[i][i];
	}

	double** backSub = Matrix::CreateMatrixDouble(size, 1);

	for (int i = size - 1; i >= 0; i--)
	{
		double sum = 0;
		for (int j = size - 1; j > i; j--)
		{
			sum += choleskyTransp.GetMatrixDouble()[i][j] * backSub[j][0];
		}
		backSub[i][0] = (forwardSub[i][0] - sum) / choleskyTransp.GetMatrixDouble()[i][i];
	}

	Matrix::DestroyMatrixDouble(forwardSub, size);

	return Matrix(backSub, size, 1);
}

//<summary>Adds a matrix to the bottom of another matrix </summary>
//<ori>Original matrix </ori>
//<toAdd>Matrix to be added to the original matrix </toAdd>
//<comment>Both matrices must have the same number of columns. </comment>
Matrix MatrixOperation::AddMatrixBottom(Matrix &ori, Matrix &toAdd) {
	int newSizeX = ori.GetDimX() + toAdd.GetDimX();
	int sizeY = ori.GetDimY();
	double** matrix = Matrix::CreateMatrixDouble(newSizeX, sizeY);
	for (int i = 0; i < newSizeX; i++) {
		for (int j = 0; j < sizeY; j++) {
			if (i < ori.GetDimX()) {
				matrix[i][j] = ori.GetMatrixDouble()[i][j];
			}
			else {
				matrix[i][j] = toAdd.GetMatrixDouble()[i - ori.GetDimX()][j];
			}
		}
	}
	return Matrix(matrix, newSizeX, sizeY);
}

//<summary>Adds a matrix to the right of another matrix </summary>
//<ori>Original matrix </ori>
//<toAdd>Matrix to be added to the original matrix </toAdd>
//<comment>Both matrices must have the same number of rows. </comment>
Matrix MatrixOperation::AddMatrixRight(Matrix &ori, Matrix &toAdd) {
	int newSizeY = ori.GetDimY() + toAdd.GetDimY();
	int sizeX = ori.GetDimX();
	double** matrix = Matrix::CreateMatrixDouble(sizeX, newSizeY);
	for (int i = 0; i < sizeX; i++) {
		for (int j = 0; j < newSizeY; j++) {
			if (j < ori.GetDimY()) {
				matrix[i][j] = ori.GetMatrixDouble()[i][j];
			}
			else {
				matrix[i][j] = toAdd.GetMatrixDouble()[i][j - ori.GetDimY()];
			}
		}
	}
	return Matrix(matrix, sizeX, newSizeY);
}

//<summary>Adds a matrix to a specific position of another matrix </summary>
//<ori>Original matrix </ori>
//<toAdd>Matrix to be added to the original matrix </toAdd>
//<initRowPos>Row position that the added matrix will be inserted in the original matrix, in real positions</initRowPos>
//<initColPos>Column position that the added matrix will be inserted in the original matrix, in real positions</initColPos>
//<comment>The original matrix' dimensions are not changed, so the added matrix must not surpass the original matrix' dimensions </comment>
void MatrixOperation::AddMatrixAtPosition(Matrix &ori, Matrix &toAdd, int initRowPos, int initColPos) {
	int nX = toAdd.GetDimX();
	int nY = toAdd.GetDimY();

	for (int i = 0; i < nX; i++) {
		for (int j = 0; j < nY; j++) {
			ori.GetMatrixDouble()[initRowPos + i - 1][initColPos + j - 1] = toAdd.GetMatrixDouble()[i][j];
		}
	}
}

//<summary>Returns a subset of a given matrix, based on the specified location </summary>
//<dim>The position, after which, all the elements will be returned as a new matrix </dim>
//<comment>This function returns a matrix with the same number of rows and the remaining number of columns from the original matrix minus the specified position. Ex. 20 - 17 = 3, the last three columns. </comment>
Matrix MatrixOperation::ExtractMatrixFromEnd(Matrix &m, int dim) {
	double** matrix = Matrix::CreateMatrixDouble(m.GetDimX(), m.GetDimY() - dim);
	for (int i = 0; i < m.GetDimX(); i++)
	{
		for (int j = m.GetDimY() - dim; j < m.GetDimY(); j++)
		{
			matrix[i][j - (m.GetDimY() - dim)] = m.GetMatrixDouble()[i][j];
		}
	}
	return Matrix(matrix, m.GetDimX(), m.GetDimY() - dim);
}

//<summary>Converts a given matrix in its normalized upper triangular form </summary>
Matrix MatrixOperation::GetNormalizedUpperTriangular(Matrix &m) {
	double** matrix = CopyMatrixDouble(m);

	for (int k = 0; k < m.GetDimX(); k++) // This index keeps zeroeing everything related to this row.
	{
		if (k < m.GetDimX() - 1)
		{
			for (int i = 0; i < m.GetDimX(); i++)
			{
				if (i == k)
				{
					double n1 = matrix[i][i];
					for (int j = 0; j < m.GetDimY(); j++)
					{
						matrix[i][j] = matrix[i][j] / n1;
					}
				}
				else if (i > k)
				{
					double n1 = matrix[i][k];
					for (int j = 0; j < m.GetDimY(); j++)
					{
						matrix[i][j] = matrix[i][j] - n1 * matrix[k][j];
					}
				}
				else
				{
					//do nothing if we are analyzing a row that has already been normalized
				}
			}
		}
		else // i.e., for the last row just divide the last element by itself to make it 1.
		{
			double n1 = matrix[k][k];
			for (int j = 0; j < m.GetDimY(); j++)
			{
				matrix[k][j] = matrix[k][j] / n1;
			}

		}
	}
	return Matrix(matrix, m.GetDimX(), m.GetDimY());
}

//<summary>Converts a given matrix in its normalized lower triangular form </summary>
Matrix MatrixOperation::GetNormalizedLowerTriangular(Matrix &m) {
	double** matrix = CopyMatrixDouble(m);

	for (int k = m.GetDimX() - 1; k >= 0; k--) // This index keeps zeroeing everything related to this row.
	{
		if (k > 0)
		{
			for (int i = m.GetDimX() - 1; i >= 0; i--)
			{
				if (i == k)
				{
					double n1 = matrix[i][i];
					for (int j = 0; j < m.GetDimY(); j++)
					{
						matrix[i][j] = matrix[i][j] / n1;
					}
				}
				else if (i < k)
				{
					double n1 = matrix[i][k];
					for (int j = 0; j < m.GetDimY(); j++)
					{
						matrix[i][j] = matrix[i][j] - n1 * matrix[k][j];
					}
				}
				else
				{
					//do nothing if we are analyzing a row that has already been normalized
				}
			}
		}
		else // i.e., for the last row just divide the last element by itself to make it 1.
		{
			double n1 = matrix[k][k];
			for (int j = 0; j < m.GetDimY(); j++)
			{
				matrix[k][j] = matrix[k][j] / n1;
			}

		}
	}
	return Matrix(matrix, m.GetDimX(), m.GetDimY());
}

//<summary>Creates a new matrix out of a set of specified columns from another matrix </summary>
//<ori>The original matrix </ori>
//<nCols>The number of columns from the original matrix that will be extracted </nCols>
//<initCol>The position of the first column to be extracted </initCol>
//<step>The space between columns that are going to be extracted</step>
//<comment>This function does not changes the original matrix</comment>
Matrix MatrixOperation::ExtractSetOfColumns(Matrix &ori, int nCols, int initCol, int step) {
	double** matrix = Matrix::CreateMatrixDouble(ori.GetDimX(), nCols);
	for (int i = 0; i < ori.GetDimX(); i++) { //for each row
		for (int j = 0; j < nCols; j++) {
			matrix[i][j] = ori.GetMatrixDouble()[i][initCol + step * j];
		}
	}
	return Matrix(matrix, ori.GetDimX(), nCols);
}

//<summary>Adds a new line and a new column with a specified term in the diagonal of this new set of rows and columns</summary>
//<ori>The original matrix </ori>
//<posLine>The position where the line will be added </posLine>
//<posCol>The position where the column will be added </posCol>
//<val>The value of the term that will be added in the crossing of the new line and new column</val>
//<comment>This function destroys the original matrix and returns the new one with the added line and column.</comment>
Matrix MatrixOperation::AddLineAndColumnWithTermAtPosition(Matrix &ori, int posLine, int posCol, double val) {
	int dimX = ori.GetDimX() + 1;
	int dimY = ori.GetDimY() + 1;
	double** matrix = Matrix::CreateMatrixDouble(dimX, dimY);
	for (int i = 0; i < dimX; i++) { //for each row
		for (int j = 0; j < dimY; j++) {
			if (i < posLine - 1)
			{
				if (j < posLine - 1) {
					matrix[i][j] = ori.GetMatrixDouble()[i][j];
				}
				else if (j == posLine - 1) {
					matrix[i][j] = 0;
				}
				else {
					matrix[i][j] = ori.GetMatrixDouble()[i][j - 1];
				}
			}
			else if (i == posLine - 1)
			{
				if (j == posLine - 1) {
					matrix[i][j] = val;
				}
				else {
					matrix[i][j] = 0;
				}
			}
			else {
				if (j < posLine - 1) {
					matrix[i][j] = ori.GetMatrixDouble()[i - 1][j];
				}
				else if (j == posLine - 1) {
					matrix[i][j] = 0;
				}
				else {
					matrix[i][j] = ori.GetMatrixDouble()[i - 1][j - 1];
				}
			}
		}
	}

	return Matrix(matrix, ori.GetDimX() + 1, ori.GetDimY() + 1);
}

//<summary>Uses the specified support and load vectors to obtain the reduced version of a complete global stiffness matrix</summary>
Matrix MatrixOperation::GetReducedMatrix(Matrix &GlobalMatrix, std::vector<Support> &vecSup, std::vector<Node> &vecNode) {
	int DOF = 6; //for shell elements
	int remove = 0; // Keeps track of the amount of support conditions to calculate the reduced size of the reduced matrix
	Matrix reduced(CopyMatrixDouble(GlobalMatrix), GlobalMatrix.GetDimX(), GlobalMatrix.GetDimY());
	for (int i = 0; i < vecSup.size(); i++) { //for each support
		int nodeID = vecSup[i].GetNode();
		for (int j = 0; j < vecSup[i].GetSupportVector().size(); j++) { //for each direction with support on the support element
			double dir = vecSup[i].GetSupportVector()[j][0]; //gets the direction
			reduced = DeleteRowAndColumn(reduced, (nodeID - 1) * DOF + (dir - 1) - remove); //deletes the line and column associated with the DOF
			remove++;
		}
	}
	return reduced;
}

void MatrixOperation::SwapLine(Matrix &matrix, int origIndex, int finalIndex) {

	for (int i = 0; i < matrix.GetDimY(); i++) //for each col
	{
		double origFinalLineTerm = matrix.GetMatrixDouble()[finalIndex][i];
		matrix.GetMatrixDouble()[finalIndex][i] = matrix.GetMatrixDouble()[origIndex][i];
		matrix.GetMatrixDouble()[origIndex][i] = origFinalLineTerm;
	}
}

void MatrixOperation::MoveLineToEnd(Matrix &matrix, int origIndex) {

	std::vector<double> origLineTerms;
	origLineTerms.reserve(matrix.GetDimY());

	for (int i = 0; i < matrix.GetDimY(); i++) //for each col, copy the original line
	{
		origLineTerms.emplace_back(matrix.GetMatrixDouble()[origIndex][i]);
	}

	for (int i = origIndex; i < matrix.GetDimX() - 1; i++) {
		for (int j = 0; j < matrix.GetDimY(); j++) {
			matrix.GetMatrixDouble()[i][j] = matrix.GetMatrixDouble()[i + 1][j];
		}
	}

	for (int i = 0; i < matrix.GetDimY(); i++) //for each col, copy the original line
	{
		matrix.GetMatrixDouble()[matrix.GetDimX() - 1][i] = origLineTerms[i];
	}
}



void MatrixOperation::MoveColumnToEnd(Matrix &matrix, int origIndex) {

	std::vector<double> origColTerms;
	origColTerms.reserve(matrix.GetDimX());

	for (int i = 0; i < matrix.GetDimX(); i++) //for each line, copy the original Column
	{
		origColTerms.emplace_back(matrix.GetMatrixDouble()[i][origIndex]);
	}

	for (int i = 0; i < matrix.GetDimX(); i++) {
		for (int j = origIndex; j < matrix.GetDimY() - 1; j++) {
			matrix.GetMatrixDouble()[i][j] = matrix.GetMatrixDouble()[i][j + 1];
		}
	}

	for (int i = 0; i < matrix.GetDimX(); i++) //for each line
	{
		matrix.GetMatrixDouble()[i][matrix.GetDimY() - 1] = origColTerms[i];
	}
}

void MatrixOperation::SwapColumn(Matrix &matrix, int origIndex, int finalIndex) {

	for (int i = 0; i < matrix.GetDimX(); i++) //for each row
	{
		double origFinalColTerm = matrix.GetMatrixDouble()[i][finalIndex];
		matrix.GetMatrixDouble()[i][finalIndex] = matrix.GetMatrixDouble()[i][origIndex];
		matrix.GetMatrixDouble()[i][origIndex] = origFinalColTerm;
	}
}

void MatrixOperation::MoveLineAndColumnToEnd(Matrix &matrix, int origIndex) {
	MoveLineToEnd(matrix, origIndex);
	MoveColumnToEnd(matrix, origIndex);
}

Matrix MatrixOperation::ExtractMatrixBasedOnLineAndColumns(Matrix &matrix, int iniLine, int nLine, int iniCol, int nCol) {
	Matrix m(nLine, nCol);

	for (int i = 0; i < nLine; i++) {
		for (int j = 0; j < nCol; j++) {
			m.GetMatrixDouble()[i][j] = matrix.GetMatrixDouble()[iniLine + i][iniCol + j];
		}
	}
	return m;
}

double MatrixOperation::GetBiggestDiagTerm(Matrix& m) {
	int dim = m.GetDimX();
	double val = 0;

	for (int i = 0; i < dim; i++) {
		if (m.GetMatrixDouble()[i][i] > val) {
			val = m.GetMatrixDouble()[i][i];
		}
	}
	return val;
}

void MatrixOperation::ConvertToEigenMatrix(Matrix &m, Eigen::MatrixXd& eigenMatrix) {

	eigenMatrix.resize(m.GetDimX(), m.GetDimY());

	for (int i = 0; i < m.GetDimX(); i++) {
		for (int j = 0; j < m.GetDimY(); j++) {
			eigenMatrix(i, j) = m.GetMatrixDouble()[i][j];
		}
	}

}

void MatrixOperation::PopulateDiagonalOnly(std::vector<double>& terms, Matrix &m) {
	for (int i = 0; i < m.GetDimX(); i++) {
		m.GetMatrixDouble()[i][i] = terms[i];
	}
}

//toFill and filler must have the same dimensions
Matrix MatrixOperation::FillMatrixBasedOnOtherMatrix(Matrix &toFill, Matrix& filler) {
	Matrix m(toFill.GetDimX(), toFill.GetDimY());

	for (int i = 0; i < m.GetDimX(); i++) {
		for (int j = 0; j < m.GetDimY(); j++) {
			if (toFill.GetMatrixDouble()[i][j] == 0) {
				m.GetMatrixDouble()[i][j] = filler.GetMatrixDouble()[i][j];
			}
			else {
				m.GetMatrixDouble()[i][j] = toFill.GetMatrixDouble()[i][j];
			}
		}
	}
	return m;
}

Matrix MatrixOperation::Sqrt(Matrix &m) {

	Eigen::MatrixXd eigenMatrix;
	ConvertToEigenMatrix(m, eigenMatrix);

	std::vector<double> srtDiag;
	srtDiag.reserve(m.GetDimX());

	Eigen::EigenSolver<Eigen::MatrixXd> s(eigenMatrix);
	Eigen::VectorXcd eigenMatrix1 = s.eigenvalues().matrix();

	Eigen::MatrixXcd eigenMatrix2 = s.eigenvectors().matrix();

	for (int i = 0; i < s.eigenvalues().size(); i++) {
		srtDiag.emplace_back(sqrt(eigenMatrix1(i).real()));
	}
	Matrix eigenValueDiag(s.eigenvalues().size(), s.eigenvalues().size());
	PopulateDiagonalOnly(srtDiag, eigenValueDiag);

	Matrix eigenVectorMatrix(s.eigenvalues().size(), s.eigenvalues().size());
	for (int i = 0; i < s.eigenvalues().size(); i++) { //for the number of eigenvalues we have
		for (int j = 0; j < s.eigenvalues().size(); j++) { //for each term of the eigenvector
			eigenVectorMatrix.GetMatrixDouble()[j][i] = eigenMatrix2(j, i).real();
		}
	}

	Matrix invEigenVector = GetInverse(eigenVectorMatrix);

	Matrix mult2 = eigenVectorMatrix * eigenValueDiag * invEigenVector;

	return mult2;
}


