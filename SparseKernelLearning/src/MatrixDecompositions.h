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
// ------------------------------ LU Storage Class ------------------------------ //
// ------------------------------------------------------------------------------ //

class LU {
	Matrix L; //Lower triangular
	Matrix U; // Upper triangular

public: 
	LU(Matrix Lmat, Matrix Umat) {
		L = Lmat;
		U = Umat;
	}

	Matrix getL() { return L; }
	Matrix getU() { return U; }
};


LU luDecomp(Matrix& mat);



// ------------------------------------------------------------------------------ //
// ------------------------------ Helper Functions ------------------------------ //
// ------------------------------------------------------------------------------ //

inline int sign(double x);