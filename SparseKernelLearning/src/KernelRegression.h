#pragma once

#include "Vector.h"
#include "Matrix.h"
#include "Kernel.h"
#include "LinearRegression.h"



class KenrnelRegressionSolutions {
	Vector weights;
	Vector sampleCoefficients;
	double intercept;

public:
	// Consturctor 
	KenrnelRegressionSolutions() {}

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
	WeightedKernel *kernel;
	KenrnelRegressionSolutions solutions;
public:
	KernelRegression(Matrix mat, Vector vec, WeightedKernel *kern) {
		trainSamples = mat;
		trainResponse = vec;
		kernel = kern;
	}

	// getter functions
	KenrnelRegressionSolutions getSolution() { return solutions; }

	// train methods 
	void train(double ridgePen);
	void train(double ridgePen, double sparsityPen);

	// predict methods
	double predict(Vector& vec);
	Vector predict(Matrix& mat);

	// Functionality for solving sparse kernel regression
	Matrix makeTMat();
	Matrix makeQMat(Matrix& tMat);
	Vector makeBetaVec(Matrix& tMat, Matrix& qMat, double ridgePen);
	
	// Convergence functionality 
	double meanSquaredError();

	double objectiveVal(double ridgePen, double sparsityPen);
	//double objectiveVal(Matrix& centeredKernelMat);

};



// ------------------------------------------------------------------------------ //
// ------------------------------ Helper fucntions ------------------------------ //
// ------------------------------------------------------------------------------ //

double quadFormObjectiveVal(Vector& weights, Matrix& qMat, Vector& betaVec, double sparsityPen);

Vector minimizeQuadForm(Matrix& qMat, Vector& betaVec, Vector& weights,  double sparsityPen, int maxIter = 1e5, double convgThresh = 1e-10);
