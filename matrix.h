#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <optional>
#include <utility>

class Matrix {
public:
	Matrix();
	Matrix(int row, int column);

	Matrix& operator=(const Matrix& other);
	Matrix operator*(const Matrix& other);
	Matrix operator+(const Matrix& other);
	Matrix operator-(const Matrix& other);
		
	std::pair<int, int> GetSize() const;
	const double& Num(int row, int column) const;
	double& Num(int row, int column);
	const std::vector<std::vector<double>>& GetMat() const;
	std::vector<std::vector<double>>& GetMat();
	void Fin(int r, int c);
	double GetTrace() const;
	double GetDeterminant() const;
	int GetRank();
	bool LIFileT(const std::string& directory) const;
	bool LIFileB(const std::string& directory) const;
	bool RFFileT(const std::string& directory);
	bool RFFileB(const std::string& directory);

	friend Matrix AtUnitMatrix(Matrix m);
	friend Matrix operator*(const Matrix& l, double num);

protected:
	int r_count = 0, col_count = 0;
	std::vector<std::vector<double>> mat;
	double det;
	int rank;
};

class UnitMatrix : public Matrix {
public:
	UnitMatrix(int size) : Matrix(size, size) {
		for (int i = 0; i < size; i++) {
			this->Num(i, i) = 1;
		}
	}
};

class DiagonalMatrix : public UnitMatrix {
public:
	DiagonalMatrix(int size) : UnitMatrix(size) {}
};

class UpperTriangularMatrix : public Matrix {
public:
	UpperTriangularMatrix(int size) : Matrix(size, size) {}
};

class LowerTriangularMatrix : public Matrix {
	LowerTriangularMatrix(int size) : Matrix(size, size) {}
};

class SymmetricMatrix : public Matrix {
	SymmetricMatrix(int size) : Matrix(size, size) {}
};

std::ostream& operator<<(std::ostream&, const Matrix& mat);
std::istream& operator>>(std::istream&, Matrix& mat);
bool operator==(const Matrix& l, const Matrix& r);

double ScalarP(const Matrix& l, const Matrix& r);
Matrix HadamardP(const Matrix& l, const Matrix& r);
double EuclidianNorm(const Matrix& m);
double MaxNorm(const Matrix& m);
double FrobeniusNorm(const Matrix& m);
double GetAngle(const Matrix& l, const Matrix& r);
Matrix InverseMatrix(const Matrix& m);
void Transpose(Matrix& m);