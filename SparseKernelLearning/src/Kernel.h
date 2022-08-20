#pragma once

#include "Vector.h"
#include "Matrix.h"

// ----------------------------------------------------------------------------------- //
// ------------------------------ Weighted Kernel Class ------------------------------ //
// ----------------------------------------------------------------------------------- //


class WeightedKernel {
	int dimension;
	Vector weights; 
public:

	// Constructors 
	WeightedKernel() {} 

	WeightedKernel(int dim) { 
		dimension = dim;
		weights = onesVec(dim);
	}

	WeightedKernel(Vector weightVec) {
		dimension = weightVec.getDim();
		weights = weightVec;
	}

	WeightedKernel(int dim, Vector weightVec) {
		dimension = dim;
		weights = weightVec;
	}

	// getter and setter for weights 
	Vector getWeights() { return weights; }
	void setWeights(Vector newWeights) { weights = newWeights;}

};



// ------------------------------------------------------------------------------ //
// ------------------------------ Specific Kernels ------------------------------ //
// ------------------------------------------------------------------------------ //



class GaussianKernel : public WeightedKernel {
	double scale;
public:
	GaussianKernel(int dim, double scaleParam) : WeightedKernel(dim) {
		scale = scaleParam;
	}

	// evaluation functions 
	double eval(Vector& vec1, Vector& vec2);
	
	Vector eval(Matrix& mat, Vector& vec);

	Matrix eval(Matrix& mat1, Matrix& mat2);

	Matrix eval(Matrix& mat);

};
