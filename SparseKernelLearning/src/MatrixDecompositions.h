#pragma once

#include "Matrix.h"
#include "Vector.h"

#include <iostream>
using namespace std;

// ------------------------------------------------------------------------------ //
// ------------------------------ QR Storage Class ------------------------------ //
// ------------------------------------------------------------------------------ //


class QR {
	Matrix Q;
	Matrix R;

public:
	QR(Matrix Qmat, Matrix Rmat) {
		Q = Qmat;
		R = Rmat;
	}

	Matrix getQ() { return Q; }
	Matrix getR() { return R; }
};





QR qrDecomp(Matrix mat);

Matrix Householder(Matrix &mat);






// ------------------------------------------------------------------------------ //
// ------------------------------ Helper Functions ------------------------------ //
// ------------------------------------------------------------------------------ //

inline int sign(double x);