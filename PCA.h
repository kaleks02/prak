#pragma once

#include <vector>
#include "matrix.h"

using namespace std;

struct Nipals {
	Matrix T, P, E;
	bool it;
};

class PCA {
private:
	Matrix X, X_cen, X_cal;
	Nipals nip;
public:
	PCA(const Matrix& matrix);
	void Calibration();
	void Centering();
	const Matrix& GetPCA() const {
		return X;
	}
	Matrix GetXcen() const;
	Matrix GetXcal() const;
	Nipals Nip();

	vector<double> Scope() const;
	vector<double> Deviation() const;
	vector<double> Dispersion();
	double FullDispersion();
	double ExplainedDispersion();
	double AvgDistance();
};