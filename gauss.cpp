#include "gauss.h"

using namespace std;

void FGauss(Matrix& m, int pos, int& count) {
	auto rows = m.GetSize().first;
	auto columns = m.GetSize().second;
	auto& matrix = m.GetMat();

	int sw = 1;
	double max_num;
	int max_index = -1;

	for (int j = pos; j < rows; j++) {
		if (sw == 1) {
			max_num = m.Num(j, pos);
			max_index = j;
			sw = 0;
		}
		else {
			if (m.Num(j, pos) > max_num) {
				max_num = m.Num(j, pos);
				max_index = j;
			}
		}
	}

	if (pos != max_index) {
		std::swap(matrix[pos], matrix[max_index]);
		++count;
	}
	if (matrix[pos][pos] != 0) {
		for (int j = pos + 1; j < rows; j++) {
			double coeff = matrix[j][pos] / matrix[pos][pos];
			for (int k = 0; k < columns; k++) {
				matrix[j][k] -= matrix[pos][k] * coeff;
				if (abs(matrix[j][k] - int(abs(matrix[j][k]))) < 0.00001)
					matrix[j][k] = int(matrix[j][k]);
			}
		}
	}
}

void BGauss(Matrix& m, int pos) {
	auto& matrix = m.GetMat();
	auto rows = m.GetSize().first;
	auto columns = m.GetSize().second;
	if (matrix[pos][pos] != 0) {
		for (int j = pos - 1; j >= 0; j--) {
			double coeff = matrix[j][pos] / matrix[pos][pos];
			for (int k = 0; k < columns; k++) {
				matrix[j][k] -= matrix[pos][k] * coeff;
				if (abs(matrix[j][k] - int(abs(matrix[j][k]))) < 0.00001)
					matrix[j][k] = int(matrix[j][k]);
			}
		}
	}
}

int Gauss(Matrix& m) {
	int count = 0;
	auto rows = m.GetSize().first;
	auto columns = m.GetSize().second;
	for (int i = 0; i < min(rows, columns); i++)
		FGauss(m, i, count);
	for (int i = min(rows, columns) - 1; i >= 0; i--)
		BGauss(m, i);
	return count;
}