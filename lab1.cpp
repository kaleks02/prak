#include <iostream>
#include <iomanip>
#include <string>

#include "Matrix.h"

using namespace std;

Matrix operator*(const Matrix& l, double num) {
	auto rows = l.GetSize().first;
	auto columns = l.GetSize().second;

	Matrix result(rows, columns);
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < columns; j++) {
			result.Num(i, j) = l.Num(i, j) * num;
		}
	}
	return result;
} // умножение на число

Matrix Matrix::operator+(const Matrix& r) {
	auto& l = *this;
	auto rows = l.GetSize().first;
	auto columns = l.GetSize().second;

	if (rows != r.GetSize().first || columns != r.GetSize().second) {
		throw runtime_error("Different sizes");
	}

	Matrix result(rows, columns);
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < columns; j++) {
			result.Num(i, j) = l.Num(i, j) + r.Num(i, j);
		}
	}
	return result;
} // сложение

Matrix Matrix::operator-(const Matrix& r) {
	auto& l = *this;
	return l + r * (-1.0);
} // вычитание

Matrix Matrix::operator*(const Matrix& r) {
	auto& l = *this;
	auto left_rows = l.GetSize().first;
	auto left_columns = l.GetSize().second;
	auto right_rows = r.GetSize().first;
	auto right_columns = r.GetSize().second;

	if (left_columns != right_rows) {
		throw runtime_error("wrong size");
	}

	Matrix result(left_rows, right_columns);
	for (int i = 0; i < left_rows; i++) {
		for (int j = 0; j < right_columns; j++) {
			double cmp = 0;
			for (int k = 0; k < left_columns; k++) {
				cmp += l.Num(i, k) * r.Num(k, j);
			}
			if (abs(cmp) - int(abs(cmp)) < 0.00001)
				cmp = int(cmp);
			result.Num(i, j) = cmp;
		}
	}
	return result;
} // перемножение матриц

Matrix HadamardP(const Matrix& l, const Matrix& r) {
	auto left_rows = l.GetSize().first;
	auto left_columns = l.GetSize().second;
	auto right_rows = r.GetSize().first;
	auto right_columns = r.GetSize().second;
	if (left_rows != right_rows || left_columns != right_columns) {
		throw runtime_error("Different size");
	}

	Matrix result(left_rows, left_columns);
	for (int i = 0; i < left_rows; i++) {
		for (int j = 0; j < left_columns; j++) {
			result.Num(i, j) = l.Num(i, j) * r.Num(i, j);
		}
	}

	return result;
} //ѕроизведение јдамара

istream& operator>>(istream& ic, Matrix& mat) {
	int rows, columns;
	ic >> rows >> columns;
	if (rows < 0 || columns < 0) {
		throw std::out_of_range("rows and columns must be non-negative");
	}

	if (rows == 0 || columns == 0) {
		rows = columns = 0;
	}

	mat.Fin(rows, columns);
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < columns; j++) {
			double current_num;
			ic >> current_num;
			mat.Num(i, j) = current_num;
		}
	}
	return ic;
}

ostream& operator<<(ostream& oc, const Matrix& mat) {
	auto rows = mat.GetSize().first;
	auto columns = mat.GetSize().second;
	oc << rows << (char)' ' << columns << '\n';
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < columns; j++) {
			oc << fixed << setprecision(3) << mat.Num(i, j) << "   ";
		}
		oc << '\n';
	}
	return oc;
}

bool operator==(const Matrix& l, const Matrix& r) {
	int left_rows = l.GetSize().first;
	int left_columns = l.GetSize().second;
	int right_rows = r.GetSize().first;
	int right_columns = r.GetSize().second;
	if (left_rows != right_rows || left_columns != right_columns) {
		throw false;
	}
	for (int i = 0; i < left_rows; i++) {
		for (int j = 0; j < left_columns; j++) {
			if (l.Num(i, j) != r.Num(i, j))
				return false;
		}
	}
	return true;
}

bool operator!=(const Matrix& l, const Matrix& r) {
	return !(l == r);
}

Matrix& Matrix::operator=(const Matrix& other) {
	if (this != &other) {
		this->det = 0;
		this->col_count = other.col_count;
		this->r_count = other.r_count;
		this->mat = other.mat;
	}
	return *this;
}