#pragma once

#include "matrix.h"
#include <cmath>
#include <algorithm>

void FGauss(Matrix& m, int pos, int& SwapCount);

void BGauss(Matrix& m, int pos);

int Gauss(Matrix& m);
