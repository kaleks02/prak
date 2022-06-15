#include <cmath>
#include <iomanip>
#include "matrix.h"
#include "gauss.h"

using namespace std;

bool IsNVector(const vector<double>& v) {
	int count = 0;
	for (int i = 0; i < v.size(); ++i) {
		if (v[i] == 0) {
			++count;
		}
	}
	return count == v.size();
}

double GetAngle(const Matrix& l, const Matrix& r) {
	return acos(ScalarP(l, r) / (EuclidianNorm(l) * EuclidianNorm(r)));
}

int Matrix::GetRank(){
	if (rank != 0) {
		int r = 0;
		Matrix cmp(*this);
		Gauss(cmp);
		auto& mat = cmp.GetMat();
		for (const auto& line : mat) {
			if (!IsNVector(line))
				++r;
		}
		rank = r;
	}
	return rank;
}

Matrix AtUnitMatrix(Matrix m) {
	auto rows = m.GetSize().first;
	auto columns = m.GetSize().second;
	if (rows != columns)
		throw runtime_error("A unit matrix cannot be attached to a non-square matrix");
	auto& matrix = m.GetMat();
	
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < columns; j++) {
			if (i == j)
				matrix[i].push_back(1);
			else 
				matrix[i].push_back(0);
		}
	}
	m.col_count *= 2;
	return m;
}

Matrix InverseMatrix(const Matrix& m) {
	auto rows = m.GetSize().first;
	auto columns = m.GetSize().second;

	if (rows != columns)
		throw runtime_error("Only a square matrix has an inverse matrix");
	if (m.GetDeterminant() == 0)
		throw runtime_error("The inverse matrix only happens for a matrix with a non-zero determinant");

	auto AtMatrix = AtUnitMatrix(m);
	Gauss(AtMatrix);
	for (int i = 0; i < rows; i++) {
		double coeff;
		for (int j = 0; j < columns * 2; j++) {
			if (i == j)
				coeff = AtMatrix.Num(i, j);
			if (j >= columns)
				AtMatrix.Num(i, j) /= coeff;
		}
	}
	Matrix Inverse(rows, rows);
	for (int i = 0; i < rows; i++) {
		for (int j = columns; j < columns * 2; j++) {
			Inverse.Num(i, j - columns) = AtMatrix.Num(i, j);
		}
	}
	return Inverse;
}

void Transpose(Matrix& m) {
	auto rows = m.GetSize().first;
	auto columns = m.GetSize().second;
	Matrix transp(columns, rows);
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < columns; j++) {
			transp.Num(j, i) = m.Num(i, j);
		}
	}
	m = transp;
}