#include "KernelRegression.h"
#include "MatrixDecompositions.h"
#include <iostream>
#include <cmath>

using namespace std;

// ------------------------------------------------------------------------------- //
// ------------------------------ Kernel Regression ------------------------------ //
// ------------------------------------------------------------------------------- //


KenrnelRegressionSolutions KernelRegression::train(double ridgePen) {
	int numSamples = trainSamples.getNumRows(); 

	// Center response
	double responseMean = mean(trainResponse);
	Vector transformedVec = trainResponse - responseMean;

	// Form kernel matrix and center rows and columns 
	Matrix kernelMatrix = kernel->eval(trainSamples);
	ColRowMeans kernColRowMeans = centerRowsCols(kernelMatrix); // centers kernel matrix 


	// add ridge penalty diagonal
	for (int coord = 0; coord < numSamples; coord++) {
		kernelMatrix.values[coord][coord] += numSamples * ridgePen;
	}

	Matrix reflectionVecs = Householder(kernelMatrix);

	// Transform response using reflection vectors

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

KenrnelRegressionSolutions KernelRegression::train(double ridgePen, double sparsityPen) {
	cout << "\n";
	solutions = this->train(ridgePen);

	int iterNum = 0; 
	double error = 1.0;
	double oldObjVal = 0;

	do {
		solutions = this->train(ridgePen);
		Vector viewWeights = solutions.getWeights();
		cout << viewWeights << "\n\n";

		Matrix tmat = this->makeTMat();
		Matrix Q = this->makeQMat(tmat);
		Matrix kernMat = kernel->eval(trainSamples);
		Vector beta = this->makeBetaVec(kernMat, tmat, Q, ridgePen);
		Vector weights = this->minimizeQuadForm(Q, beta, sparsityPen / 2);
		solutions.setWeights(weights);
		kernel->setWeights(weights);

		double newObjVal = this->objectiveVal(ridgePen, sparsityPen);
		error = abs(newObjVal - oldObjVal);
		oldObjVal = newObjVal;
		iterNum++;

		cout << iterNum << " " << oldObjVal << endl;

	} while ((iterNum < 100) && (error > 1e-10));

	return solutions;
}


double KernelRegression::predict(Vector& vec) {
	Vector coeffs = getSolution().getSampleCoefficients();
	double coeffMean = mean(coeffs);
	coeffs = coeffs - coeffMean;
	Vector kernVec = kernel->eval(trainSamples, vec);
	double kernDotProd = coeffs * kernVec;

	double prediction = getSolution().getIntercept() + kernDotProd;
	return prediction;
}

Vector KernelRegression::predict(Matrix& mat) {
	int nrows = mat.getNumRows();
	//Vector predictions(nrows);
	Vector coeffs = getSolution().getSampleCoefficients();
	double coeffMean = mean(coeffs);
	coeffs = coeffs - coeffMean;
	Matrix kernMat = kernel->eval(trainSamples, mat);
	Vector kernDotProds = kernMat * coeffs;
	double intercept = getSolution().getIntercept();
	return kernDotProds + intercept;
}



// Functionality for sparse kernel regression 
Matrix KernelRegression::makeTMat() {
	int nrows = trainSamples.getNumRows(), ncols = trainSamples.getNumCols();
	Matrix tmat = zeroMat(nrows, ncols);
	Matrix derivMat(nrows, ncols);
	Vector coeffs = this->getSolution().getSampleCoefficients();

	// Center coefficients 
	double coeffMean = mean(coeffs);
	coeffs = coeffs - coeffMean;

	for (int row = 0; row < nrows; row++) {
		Vector rowVec = trainSamples.getRow(row);
		derivMat = kernel->deriv(trainSamples, rowVec);

		derivMat = coeffs.values[row] * derivMat;
		tmat = tmat + derivMat;
	}
	return tmat;
}

Matrix KernelRegression::makeQMat(Matrix& tMat) {
	int ncols = tMat.getNumCols(), nrows = tMat.getNumRows();
	Matrix Q = crossprod(tMat);
	Q = Q / nrows;

	// centert rows and columns  
	ColRowMeans qMeans = centerRowsCols(Q);
	return Q; 
}

Vector KernelRegression::makeBetaVec(Matrix& kernMat,  Matrix& tMat, Matrix& qMat, double ridgePen) {
	int nrows = tMat.getNumRows(), ncols = tMat.getNumCols();
	double responseMean = mean(trainResponse);
	Vector centeredResponse = trainResponse - responseMean; 
	Vector weights = this->getSolution().getWeights();
	Vector coeffs = this->getSolution().getSampleCoefficients();
	Matrix tmatTranspose = transpose(tMat);

	// column center tmat 
	Vector colMeans = tMat.colMeans();

	for (int col = 0; col < ncols; col++) {
		for (int row = 0; row < nrows; row++) {
			tMat.values[row][col] -= colMeans.values[col];
		}
	}

	// add ridge penalty to centered kernel matrix 
	for (int coord = 0; coord < nrows; coord++) {
		kernMat.values[coord][coord] += ridgePen * nrows / 2; 
	}


	double part1Scale = 1 / nrows;
	Vector kernCoeff = kernMat * coeffs;
	Vector part1 = centeredResponse - kernCoeff;
	part1 = tmatTranspose * part1;
	part1 = part1 * part1Scale; 

	Vector part2 = qMat * weights;
	Vector beta = part1 + part2;
	return beta;
}


Vector KernelRegression::minimizeQuadForm(Matrix& qMat, Vector& betaVec, double sparsityPen, int maxIter, double convgThresh) {
	double error = 1;
	int iterNum = 0;
	int dim = qMat.getNumCols();
	Vector weights = this->getSolution().getWeights();

	double oldObjVal = 0; 

	do {
		for (int coord = 0; coord < dim; coord++) {
			Vector qRow = qMat.getRow(coord);
			double dotProd = qRow* weights;
			double resid = dotProd - qRow.values[coord] * weights.values[coord];
			resid = betaVec.values[coord] - resid;
			weights.values[coord] = softThresh(resid, sparsityPen);

			// scale by qMat diagonal entry 
			weights.values[coord] = weights.values[coord] / qMat.values[coord][coord];

			// project to +/- 1 if value is above or below 
			if (weights.values[coord] > 1) { weights.values[coord] = 1; }
			else if (weights.values[coord] < -1) { weights.values[coord] = -1; }
		}

		// update error 
		iterNum++;
		double newObjVal = quadFormObjectiveVal(weights, qMat, betaVec, sparsityPen);
		error = (newObjVal - oldObjVal) / (oldObjVal + 1e-5);
		oldObjVal = newObjVal;

	} while ((iterNum <= maxIter) && (error > convgThresh));
	return weights;
}



// Convergence Functionality
double KernelRegression::meanSquaredError() {
	int nrows = trainSamples.getNumRows();

	Vector predictions = this->predict(trainSamples);
	Vector residuals = predictions - trainResponse;

	double normVal = norm(residuals);
	return normVal * normVal / nrows;
}

double KernelRegression::objectiveVal(double ridgePen, double sparsityPen) {
	double mse = this->meanSquaredError();

	// ridge penalty
	Matrix kernelMat = kernel->eval(trainSamples);
	ColRowMeans kernColRowMeans = centerRowsCols(kernelMat); // centers kernel matrix
	Vector coeffs = this->getSolution().getSampleCoefficients();
	Vector kernelCoeffs = kernelMat * coeffs;
	double ridgeObj = kernelCoeffs * coeffs;
	ridgeObj *= ridgePen;

	// sparsity penalty 
	Vector weights = this->getSolution().getWeights();
	double sparsityNorm = 0;
	for (int coord = 0; coord < weights.getDim(); coord++) {
		sparsityNorm += abs(weights.values[coord]);
	}
	sparsityNorm *= sparsityPen;

	return mse + ridgeObj + sparsityNorm;
}


// ------------------------------------------------------------------------------ //
// ------------------------------ Helper fucntions ------------------------------ //
// ------------------------------------------------------------------------------ //

double quadFormObjectiveVal(Vector& weights, Matrix& qMat, Vector& betaVec, double sparsityPen) {
	Vector qWeight = qMat * weights;
	double weightQweight = weights * qWeight;
	weightQweight /= 2;

	double betaWeight = betaVec * weights;

	double sparsityNorm = 0;
	int dim = betaVec.getDim();
	for (int coord = 0; coord < dim; coord++) {
		sparsityNorm += abs(weights.values[coord]);
	}
	sparsityNorm *= sparsityPen;

	return weightQweight - betaWeight + sparsityNorm;
}