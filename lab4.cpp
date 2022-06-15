#include "Matrix.h"
#include <fstream>

using namespace std;

bool Matrix::LIFileB(const string& file) const {
	ofstream out_str;
	out_str.open(file, std::ifstream::binary);
	if (out_str.is_open()) {
		auto rows = GetSize().first;
		auto columns = GetSize().second;
		out_str.write((char*)&rows, sizeof rows);
		out_str.write((char*)&columns, sizeof columns);
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < columns; j++)
				out_str.write((char*)&Num(i, j), sizeof(double));
	}
	return true;
}

bool Matrix::LIFileT(const string& file) const {
	ofstream out_str;
	out_str.open(file);
	if (out_str.is_open()) {
		out_str << *this;
	}
	return true;
}


bool Matrix::RFFileT(const string& file) {
	ifstream in_str;
	in_str.open(file);
	if (in_str.is_open()) {
		in_str >> *this;
	}
	return true;
}

bool Matrix::RFFileB(const string& file) {
	ifstream in_str;
	in_str.open(file, std::ifstream::binary);
	if (in_str.is_open()) {
		in_str.read((char*)&r_count, sizeof(int));
		in_str.read((char*)&col_count, sizeof(int));
		Fin(r_count, col_count);
		for (int i = 0; i < r_count; i++)
			for (int j = 0; j < col_count; j++)
				in_str.read((char*)&Num(i, j), sizeof(double));
	}
	return true;
}