#pragma once

#include "Vector.h"
#include "Matrix.h"
#include "Kernel.h"


class KenrnelRegressionSolutions {
	Vector weights;
	Vector sampleCoefficients;
	double intercept;

public:
	// Consturctor 
	KenrnelRegressionSolutions(Vector w, Vector sc, double i) {
		weights = w;
		sampleCoefficients = sc;
		intercept = i;
	}

	// Getter and setter functions 
	Vector getWeights() { return weights; }
	Vector getSampleCoefficients() { return sampleCoefficients; }
	double getIntercept() { return intercept; }

	void setWeights(Vector& newWeights) { weights = newWeights; }
	void setSampleCoefficients(Vector& newSampleCoefficients) { sampleCoefficients = newSampleCoefficients; }
	void setIntercept(double newIntercept) { intercept = newIntercept; }
};



// ------------------------------------------------------------------------------- //
// ------------------------------ Kernel Regression ------------------------------ //
// ------------------------------------------------------------------------------- //



class KernelRegression {
	Matrix trainSamples;
	Vector trainResponse;
	WeightedKernel kernel;
public:
	KernelRegression(Matrix mat, Vector vec, WeightedKernel kern) {
		trainSamples = mat;
		trainResponse = vec;
		kernel = kern;
	}

	// train methods 
	KenrnelRegressionSolutions train(double ridgePen);
	KenrnelRegressionSolutions train(double ridgePen, double sparsityPen);
};




