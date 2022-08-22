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

	// pure virtual evaluation functions 
	virtual double eval(Vector& vec1, Vector& vec2) = 0;

	virtual Vector eval(Matrix& mat, Vector& vec) = 0;

	virtual Matrix eval(Matrix& mat1, Matrix& mat2) = 0;

	virtual Matrix eval(Matrix& mat) = 0;

	// pure virtual derivative functions
	virtual Vector deriv(Vector& vec1, Vector& vec2) = 0;
	
	virtual Matrix deriv(Matrix& mat, Vector& vec) = 0;


	// getter and setters 
	int getDim() { return dimension; }
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

	// pure virtual derivative functions
	Vector deriv(Vector& vec1, Vector& vec2);

	Matrix deriv(Matrix& mat, Vector& vec);


};
