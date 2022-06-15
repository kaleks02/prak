#include "PCA.h"
#include <cmath>

PCA::PCA(const Matrix& matrix) {
	X = matrix;
}

void PCA::Centering() {
	auto rows = X.GetSize().first;
	auto columns = X.GetSize().second;

	Matrix result(rows, columns);
	for (int i = 0; i < columns; ++i) {
		double avg = 0;
		for (int j = 0; j < rows; ++j) {
			avg += X.Num(j, i) / rows;
		}
		for (int j = 0; j < rows; ++j) {
			result.Num(j, i) = X.Num(j, i) - avg;
		}
	}
	X_cen = result;
}

void PCA::Calibration() {
	auto rows = X.GetSize().first;
	auto columns = X.GetSize().second;

	Matrix result(rows, columns);
	for (int i = 0; i < columns; ++i) {
		double s2 = 0;
		double avg = 0;
		for (int j = 0; j < rows; ++j) {
			avg += X.Num(j, i) / rows;
		}
		for (int j = 0; j < rows; ++j) {
			s2 += (X.Num(j, i) - avg) * (X.Num(j, i) - avg) / ((rows - 1) * 1.0);
		}
		for (int j = 0; j < rows; j++) {	
			result.Num(j, i) = (X.Num(j, i) - avg) / sqrt(s2);
		}
	}
	X_cal = result;
}

Matrix PCA::GetXcal() const {
	return X_cal;
}

Matrix PCA::GetXcen() const {
	return X_cen;
}

Nipals PCA::Nip() {
	Centering();
	Calibration();
	Matrix Ea = X_cal;
	auto rows = X.GetSize().first;
	auto columns = X.GetSize().second;
	Matrix T(rows, columns), P(rows, columns);
	Matrix E = X_cal;

	Matrix t(rows, 1), p, t_old;
	vector <Matrix> t_V, p_V;
	for (int j = 0; j < columns; ++j) {
		for (int i = 0; i < rows; ++i) {
			t.Num(i, 0) = X.Num(i, j);
		}
		do {
			Matrix tT = t; Transpose(tT);
			p = tT * Ea;
			for (auto& vector: p.GetMat()) {
				for (auto& element : vector) {
					element /= (tT * t).Num(0, 0);
				}
			}
			Transpose(p);
			for (auto& vector : p.GetMat()) {
				for (auto& element : vector) {
					element /= EuclidianNorm(p);
				}
			}	
			Matrix pT = p; Transpose(pT);
			t_old = t;
			t = Ea * p;
			for (auto& vector : t.GetMat()) {
				for (auto& element : vector) {
					element /= (pT * p).Num(0, 0);
				}
			}
		} while (EuclidianNorm(t_old - t) > 0.0001);
		Matrix pT = p; Transpose(pT);
		Ea = Ea - t * pT;
		p_V.push_back(p);
		t_V.push_back(t);
	}

	Matrix T1(rows, t_V.size()), P1(p.GetSize().first, p_V.size());
	for (int i = 0; i < rows; ++i) {
		for (int r = 0; r < t_V.size(); ++r) {
			T1.Num(i, r) = t_V[r].Num(i, 0);
			if (i < p_V[0].GetSize().first)
				P1.Num(i, r) = p_V[r].Num(i, 0);
		}
	}
	nip.T = T1; 
	nip.P = P1; 
	nip.E = Ea;
	nip.it = true;
	return nip;
}

vector<double> PCA::Scope() const {
	if (!nip.it) throw runtime_error("Для вычисления размаха следует применить NIPALS.");
	vector<double> scope;
	for (int i = 0; i < nip.T.GetSize().first; ++i) {
		Matrix t(1, nip.T.GetSize().second);
		for (int j = 0; j < nip.T.GetSize().second; ++j) {
			t.Num(0, j) = nip.T.Num(i, j);
		}
		Matrix tT = nip.T; Transpose(tT);
		Matrix Inverse_T = InverseMatrix(tT * nip.T);
		auto h = (t * Inverse_T) * tT;
		scope.push_back(h.Num(0, 0));
	}
	return scope;
}

vector<double> PCA::Deviation() const {
	if (!nip.it) throw runtime_error("Для вычисления отклонения следует применить NIPALS.");
	vector<double> v_V;
	for (int i = 0; i < nip.E.GetSize().first; ++i) {
		double v_i = 0;
		for (int j = 0; j < nip.E.GetSize().second; ++j) {
			v_i += nip.E.Num(i, j) * nip.E.Num(i, j);
		}
		v_V.push_back(v_i);
	}
	return v_V;
}

vector<double> PCA::Dispersion() {
	if (!nip.it) throw runtime_error("Для вычисления дисперсии следует применить NIPALS.");
	vector<double> devs = Deviation();
	double cmp = nip.E.GetSize().second;
	for (auto& i : devs) {
		i /= cmp;
	}
	return devs;
}

double PCA::FullDispersion() {
	if (!nip.it) throw runtime_error("Для вычисления дисперсии следует применить NIPALS.");
	vector<double> disp = Dispersion();
	double fullDispersion = 0;
	for (auto i : disp) {
		fullDispersion += i / nip.E.GetSize().first;
	}
	return fullDispersion;
}

double PCA::ExplainedDispersion() {
	if (!nip.it) throw runtime_error("Для вычисления дисперсии следует применить NIPALS.");
	double Exp = 0;
	for (size_t i = 0; i < X_cal.GetSize().first; ++i) {
		for (size_t j = 0; j < X_cal.GetSize().second; ++j) {
			Exp += X_cal.Num(i, j) * X_cal.Num(i, j);
		}
	}
	return (1 - nip.E.GetSize().first * AvgDistance() / Exp);
}

double PCA::AvgDistance() {
	if (!nip.it) throw runtime_error("Для вычисления следует применить NIPALS.");
	vector<double> devs = Deviation();
	double avgD = 0;
	for (auto element : devs) {
		avgD += element / nip.E.GetSize().first;
	}
	return avgD;
}
