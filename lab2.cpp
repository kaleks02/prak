#include <iostream>
#include <algorithm>
#include <cmath>	
#include <iomanip>
#include "matrix.h"
#include "gauss.h"

using namespace std;

double Matrix::GetTrace() const {
	if (r_count != col_count) {
		throw runtime_error("trace of the matrix is considered only for a square matrix");
	}
	double trace = 0;
	for (int i = 0; i < r_count; i++) {
		trace += this->Num(i, i);
	}
	return trace;
}

double Matrix::GetDeterminant() const {
	if (r_count == 0 || col_count == 0)
		throw runtime_error("determinant of the matrix is considered only  for a non-empty square matrix");
	if (r_count != col_count)
		throw runtime_error("determinant of the matrix is considered only for a square matrix");
	Matrix ñmp(*this);
	double deter = 1;
	for (int i = 0; i < r_count; i++) {
		deter *= ñmp.Num(i, i);
	}
	auto det = pow(-1, Gauss(ñmp)) * deter;
	return det;
}

bool isVector(const Matrix& m) {
	return m.GetSize().first == 1 || m.GetSize().second == 1;
}

double ScalarP(const Matrix& l, const Matrix& r) {
	if (isVector(l) && isVector(r)) {
		double result = 0;
		auto left_rows = l.GetSize().first;
		auto left_columns = l.GetSize().second;
		auto right_rows = r.GetSize().first;
		auto right_columns = r.GetSize().second;
		if (max(left_rows, left_columns) != max(right_rows, right_columns))
			throw invalid_argument("the dimensions of the vectors are different");
		if (right_rows == right_columns) {
			return l.Num(0, 0) * r.Num(0, 0);
		}
		else if (left_rows > left_columns && right_rows > right_columns) {
			for (int i = 0; i < left_rows; i++)
				result += l.Num(i, 0) * r.Num(i, 0);
		}
		else if (left_rows > left_columns && right_columns > right_rows) {
			for (int i = 0; i < left_rows; i++)
				result += l.Num(i, 0) * r.Num(0, i);
		}
		else if (left_columns > left_rows && right_rows > right_columns) {
			for (int i = 0; i < left_columns; i++)
				result += l.Num(0, i) * r.Num(i, 0);
		}
		else if (left_columns > left_rows && right_columns > right_rows) {
			for (int i = 0; i < left_columns; i++)
				result += l.Num(0, i) * r.Num(0, i);
		}
		return result;
	}
	else throw runtime_error("You cannot multiply matrices");
}

double EuclidianNorm(const Matrix& m) {
	if (!isVector(m))
		throw runtime_error("not a vector is entered");
	else return sqrt(ScalarP(m, m));
}

double MaxNorm(const Matrix& m) {
	if (!isVector(m))
		throw runtime_error("not a vector is entered");
	auto rows = m.GetSize().first;
	auto columns = m.GetSize().second;
	optional<double> cur_max;
	if (rows == 1) {
		for (int i = 0; i < columns; i++) {
			if (!cur_max)
				cur_max = m.Num(0, i);
			else 
				cur_max = max(m.Num(0, i), cur_max.value());
		}
	}
	else {
		for (int i = 0; i < rows; i++) {
			if (!cur_max)
				cur_max = m.Num(i, 0);
			else 
				cur_max = max(m.Num(i, 0), cur_max.value());
		}
	}
	return cur_max.value();
}

double FrobeniusNorm(const Matrix& m) {
	double result = 0;
	auto rows = m.GetSize().first;
	auto columns = m.GetSize().second;
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < columns; ++j) {
			result += m.Num(i, j);
		}
	}
	return sqrt(result);
}