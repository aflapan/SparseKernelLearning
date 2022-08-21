#include "KernelRegression.h"
#include "MatrixDecompositions.h"
#include <iostream>

using namespace std;

// ------------------------------------------------------------------------------- //
// ------------------------------ Kernel Regression ------------------------------ //
// ------------------------------------------------------------------------------- //


KenrnelRegressionSolutions KernelRegression::train(double ridgePen) {
	int numSamples = trainSamples.getNumRows(); 

	// Center response
	double responseMean = mean(trainResponse);
	trainResponse = trainResponse - responseMean;

	// Form kernel matrix and center rows and columns 
	Matrix kernelMatrix = kernel->eval(trainSamples);
	ColRowMeans kernColRowMeans = centerRowsCols(kernelMatrix); // centers kernel matrix 


	// add ridge penalty diagonal
	for (int coord = 0; coord < numSamples; coord++) {
		kernelMatrix.values[coord][coord] += numSamples * ridgePen;
	}

	Matrix reflectionVecs = Householder(kernelMatrix);

	// Transform response using reflection vectors

	Vector transformedVec = trainResponse;

	for (int index = 0; index < numSamples; index++) {
		Vector reflectVec = reflectionVecs.getCol(index);
		double projection = (reflectVec * transformedVec);
		reflectVec = (2 * projection) * reflectVec;
		transformedVec = transformedVec - reflectVec;
	}
	
	// Solve for coefficients and intercept 
	Vector coefficients = backSolve(kernelMatrix, transformedVec);
	double coeffMean = mean(coefficients);
	Vector centeredCoeffs = coefficients - coeffMean;

	Vector colMeans = kernColRowMeans.getColMeans();
	double intercept = responseMean - (centeredCoeffs* colMeans);

	KenrnelRegressionSolutions trainedSolutions(kernel->getWeights(), coefficients, intercept);

	solutions = trainedSolutions;
	return solutions;
}


double KernelRegression::predict(Vector& newSample) {
	Vector coeffs = getSolution().getSampleCoefficients();
	double coeffMean = mean(coeffs);
	coeffs = coeffs - coeffMean;
	Vector kernVec = kernel->eval(trainSamples, newSample);
	double kernDotProd = coeffs * kernVec;

	double prediction = getSolution().getIntercept() + kernDotProd;
	return prediction;
}




// ------------------------------------------------------------------------------ //
// ------------------------------ Helper fucntions ------------------------------ //
// ------------------------------------------------------------------------------ //



