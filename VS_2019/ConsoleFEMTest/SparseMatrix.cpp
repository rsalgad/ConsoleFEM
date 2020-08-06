#include "pch.h"
#include "SparseMatrix.h"
#include <map>

SparseMatrix::SparseMatrix()
{
}

SparseMatrix::SparseMatrix(const int dimX, const int dimY)
{
	_dimX = dimX;
	_dimY = dimY;
}

SparseMatrix::~SparseMatrix()
{
}

std::map<int, std::map<int, double>>* SparseMatrix::GetMap()
{
	return &_matrix;
}

void SparseMatrix::SetSize(const int size)
{
	_size = size;
}

int SparseMatrix::SparseSize(const Matrix& m)
{
	double countNonZero = 0;
	for (int i = 0; i < m.GetDimX(); i++) {
		for (int j = 0; j < m.GetDimY(); j++) {
			if (m.GetMatrixDouble()[i][j] != 0) {
				countNonZero++;
			}
		}
	}
	return countNonZero;
}

double** SparseMatrix::CreateMatrixDouble(const int size) {
	double** m = new double* [3];
	for (int i = 0; i < 3; i++)
		m[i] = new double[size] {}; // '{}' initializes everything to zero
	return m;
}

SparseMatrix SparseMatrix::ConvertToSparseMatrix(const Matrix& m)
{
	SparseMatrix matrix(m.GetDimX(), m.GetDimY());
	std::vector<double> rows, cols, vals;

	for (int i = 0; i < m.GetDimX(); i++) {
		for (int j = 0; j < m.GetDimY(); j++)
		{
			if (m.GetMatrixDouble()[i][j] != 0) { //This should be ordered by rows then cols.
				rows.emplace_back((double)i);
				cols.emplace_back((double)j);
				vals.emplace_back(m.GetMatrixDouble()[i][j]);
			}
		}
	}

	double** d = new double* [rows.size()];
	for (int i = 0; i < rows.size(); i++)
		d[i] = new double[3]{ rows[i], cols[i], vals[i] };

	//matrix.SetMatrixDouble(d);
	matrix.SetSize(rows.size());

	return matrix;
}

SparseMatrix SparseMatrix::ConvertToSymmetricSparseMatrix(const Matrix& m)
{
	SparseMatrix matrix(m.GetDimX(), m.GetDimY());
	int count = 0;

	for (int i = 0; i < m.GetDimX(); i++) {
		std::map<int, double> mapRow;
		for (int j = 0; j <= i; j++)
		{
			if (m.GetMatrixDouble()[i][j] != 0) { //This should be ordered by rows then cols.
				mapRow.insert(std::pair<int, double>(j, m.GetMatrixDouble()[i][j]));
				count++;
			}
		}
		matrix.GetMap()->insert(std::pair<int, std::map<int,double>>(i, mapRow));
	}

	matrix.SetSize(count);

	return matrix;
}

int SparseMatrix::IndexOfFirstLine(const int line) {
	for (int i = 0; i < _size; i++)
	{
		if (_matrix[i][0] == line)
			return i;
	}
	return -1;
}

int SparseMatrix::ElementsInLineXBeforeColY(const int lineIndex, const int col) {
	int index = lineIndex;
	while (_matrix[index][1] != col) {
		index++;
	}
	return index - lineIndex;
}
