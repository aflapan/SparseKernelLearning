#pragma once

#include <iostream>
#include <cmath>
#include "Vector.h"
#include "Matrix.h"
using namespace std;


// ------------------------------------------------------------------------------- //
// ------------------------------ Base Kernel Class ------------------------------ //
// ------------------------------------------------------------------------------- //

class Kernel {
	int dimension;
public:
	// default constructor
	Kernel() {};
	~Kernel() {};

	// virtual Kernel evaltuation functions. 
	virtual double eval(Vector& vec1, Vector& vec2);
	virtual Vector eval(Matrix& mat, Vecto& vec);
	virtual Matrix eval(Matrix& mat1, Matrix& mat2);
	virtual Matrix eval(Matrix& mat);
};

class WeightedKernel : public Kernel {
	Vector weights; 
public:
	// constructord
	WeightedKernel() : Kernel() { weights = onesVec(dimension); }
	WeightedKernel(Vector w) : Kernel() { weights = w; }

};



