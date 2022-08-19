#pragma once

#include "Vector.h"
#include "Matrix.h"

// ------------------------------------------------------------------------------------------------------------ //
// ------------------------------ Kernel and Weighted Kernel Polymorphic Classes ------------------------------ //
// ------------------------------------------------------------------------------------------------------------ //

class Kernel { // note to self, think about what I want to encapsualte into kernel 

public:

};


class WeightedKernel :public Kernel {
	int dimension;
	Vector weights; 
public:

	// Constructors 
	WeightedKernel(int dim) : Kernel() { 
		dimension = dim;
		weights = onesVec(dim);
	}

	WeightedKernel(Vector weightVec) : Kernel() {
		dimension = weightVec.getDim();
		weights = weightVec;
	}

	WeightedKernel(int dim, Vector weightVec) : Kernel() {
		dimension = dim;
		weights = weightVec;

	}

	// getter and setter for weights 
	Vector getWeights() { return weights; }
	void setWeights(Vector newWeights) { weights = newWeights;}

	// evaluation functions 




	// Derivative functions with respect to weights 



};



// ------------------------------------------------------------------------------ //
// ------------------------------ Specific Kernels ------------------------------ //
// ------------------------------------------------------------------------------ //



class GaussianKernel : public Kernel {
	double scale;
public:
	GaussianKernel(double scaleParam) : Kernel() {
		scale = scaleParam;
	}

	// evaluation functions 
	double eval(Vector& vec1, Vector& vec2);
	
	Vector eval(Matrix& mat, Vector& vec);

	Matrix eval(Matrix& mat1, Matrix& mat2);

	Matrix eval(Matrix& mat);

};
