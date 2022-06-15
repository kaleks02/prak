#include "matrix.h"

using namespace std;

Matrix::Matrix() : r_count(0), col_count(0) {}

Matrix::Matrix(int row, int column) : r_count(row), col_count(column) {
	if (row < 0 || column < 0) {
		throw std::out_of_range("rows and columns must be non-negative");
	} 

	if (row == 0 || column == 0) {
		row = column = 0;
	}

	r_count = row;
	col_count = column;
	mat.assign(row, std::vector<double>(column, 0));
}

std::pair<int, int> Matrix::GetSize() const {
	return std::make_pair(r_count, col_count);
}

const double& Matrix::Num(int row, int column) const {
	return mat[row][column];
}

double& Matrix::Num(int row, int column) {
	return mat[row][column];
}

const std::vector<std::vector<double>>& Matrix::GetMat() const {
	return mat;
}

std::vector<std::vector<double>>& Matrix::GetMat() {
	return mat;
}

void Matrix::Fin(int r, int c) {
	if (r < 0 || c < 0)  {
		throw std::out_of_range("rows and columns must be non-negative");
	}

	if (r == 0 || c == 0) {
		r = c = 0;
	}

	r_count = r;
	col_count = c;
	mat.assign(r, std::vector<double>(c, 0));
}