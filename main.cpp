#include <iostream>
#include <cmath>
#include "matrix.h"
#include "PCA.h"
#include <cassert>

using namespace std;
int main() {
	/*Matrix m(2, 2);
	m.Num(0, 0) = 2;
	m.Num(0, 1) = 3;
	m.Num(1, 0) = 13;
	m.Num(1, 1) = 11;
	double t = m.GetTrace();
	cout << t << endl << endl;
	Matrix m1 = InverseMatrix(m);
	cout << m1;
	cout << m * m1;*/
	/*Matrix m(3, 2);
	m.Num(0, 0) = 2;
	m.Num(0, 1) = 3;
	m.Num(1, 0) = 13;
	m.Num(1, 1) = 11;
	m.Num(2, 0) = 5;
	m.Num(2, 1) = 3;
	Transpose(m);
	m.LIFileB("main.bin");
	Matrix m2;
	m2.RFFileB("main.bin");
	cout << m2;*/
	/*Matrix v1(2, 1);
	Matrix v2(1, 2);
	v1.Num(0, 0) = 1;
	v1.Num(1, 0) = 3;
	v2.Num(0, 0) = 5;
	v2.Num(0, 1) = 4;
	double pi = 3.1415926535;
	cout << GetAngle(v1, v2) * 180 / pi*/;
	/*Matrix m(2, 2);
	m.Num(0, 0) = 2;
	m.Num(0, 1) = 3;
	m.Num(1, 0) = 13;
	m.Num(1, 1) = 11;
	cout << m.GetDeterminant();*/
	Matrix m;
	m.RFFileT("main.txt");
	PCA x(m);
	x.Centering();
	x.Calibration();
	Nipals t = x.Nip();
	cout << t.T << endl;
	/*vector <double> v = x.Scope();
	for (auto el : v) {
		cout << el << " ";
	}*/
}