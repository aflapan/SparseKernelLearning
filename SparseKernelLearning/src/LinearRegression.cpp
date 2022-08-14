/*
A class and functionality for linear models.
The solver methods allow for both Ridge Regression and the Lasso/elastic net.
Uses the Matrix and Vector classes from VectorAndMatrixClasses.cpp
*/

#include <iostream>
#include "Vector.h"
#include "Matrix.h"
#include "MatrixDecompositions.h"
#include "LinearRegression.h"
using namespace std;



Vector backSolve(Matrix& upperTriangMat, Vector& targetVec) {
	// Include check for zero on diagonal
	int dim = targetVec.getDim();
	Vector solution(dim);

	for (int coord = dim - 1; coord >= 0; coord--) {
		double knownLinearComponent = 0;
		for (int i = coord + 1; i < dim; i++) {
			knownLinearComponent += solution.values[i] * upperTriangMat.values[coord][i];
		}
		solution.values[coord] = (targetVec.values[coord] - knownLinearComponent) / upperTriangMat.values[coord][coord];
	}
	return solution;
}


LinearSolutions LinearModel::train() {

	// Normalize columns to be mean zero and standDev 1
	Vector means = trainSamples.colMeans();
	Vector stds = trainSamples.colStanDevs();
	Normalize(trainSamples, means, stds);

	double responseMean = mean(trainResponse);
	trainResponse = trainResponse - responseMean;

	QR qr = qrDecomp(trainSamples);
	Matrix Q = qr.getQ(), R = qr.getR();

	// Figure out how to make this efficient.
	// Why can't I use expressions in place of matrix objects? 
	Matrix Qt = transpose(Q);
	Vector transposeResponse = Qt * trainResponse;
	Vector regCoefficients = backSolve(R, transposeResponse);
	
	// Put coefficients on feature scale and make intercept
	int dim = regCoefficients.getDim();

	for (int coord = 0; coord < dim; coord++) {
		if (abs(stds.values[coord]) > 1e-10)
			regCoefficients.values[coord] /= stds.values[coord];
	}

	double intercept = responseMean - regCoefficients * means;
	LinearSolutions solutions(intercept, regCoefficients);

	return solutions;
}

LinearSolutions LinearModel::train(double ridgePen, double sparsityPen) {
	int nrows = trainSamples.getNumRows(), ncols = trainSamples.getNumCols();

	// Normalize columns to be mean zero and standDev 1
	Vector means = trainSamples.colMeans();
	Vector stds = trainSamples.colStanDevs();
	Normalize(trainSamples, means, stds);

	double responseMean = mean(trainResponse);
	trainResponse = trainResponse - responseMean;

	// Initialize solutions 
	Vector coefficients = zeroVec(ncols);
	double intercept = 0;
	LinearSolutions solutions(intercept, coefficients);

	// Initialize vaiables
	double shrinkage = nrows * sparsityPen / 2;
	int iterNum = 0;
	double error = 1.0;
	double oldObjVal = this->objectiveValue(solutions, ridgePen, sparsityPen);

	// Form X^tX matrix with ridge diagonal 
	Matrix transposeData = transpose(trainSamples);
	Matrix XtX = transposeData * trainSamples;
	for (int index = 0; index < ncols; index++) {
		XtX.values[index][index] += ridgePen * nrows;
	}

	while ((error > 1e-10) && (iterNum < 1e5)) {
		// Update all coefficients 
		for (int coeffIndex = 0; coeffIndex < ncols; coeffIndex++) {
			Vector dataCol = trainSamples.getCol(coeffIndex);
			
			Vector rowVec = XtX.getRow(coeffIndex);
			double resid = rowVec * coefficients - coefficients.values[coeffIndex] * rowVec.values[coeffIndex]; 
			double updatedCoeff = softThresh(dataCol * trainResponse - resid, shrinkage);
			coefficients.values[coeffIndex] = updatedCoeff / rowVec.values[coeffIndex];
			solutions.updateCoeffVal(coefficients.values[coeffIndex], coeffIndex);
		}

		double newObjVal = this->objectiveValue(solutions, ridgePen, sparsityPen);
		error = abs(oldObjVal - newObjVal);
		iterNum++;
		oldObjVal = newObjVal;
	}

	return solutions;
}



ostream& operator<<(ostream& stream, LinearSolutions& solution) {
	Vector coeffs = solution.getCoefficients();
	int dim = coeffs.getDim();

	stream << "Intercept: " << solution.getIntercept() << "\n\n";

	stream << "Regression Coefficients: " << "\n";
	for (int coord = 0; coord < dim; coord++) {
		stream << coeffs.values[coord] << ' ';
	}
	stream << "\n";
	return stream;
}

// ------------------------------------------------------------------------------ //
// ------------------------------ Helper Functions ------------------------------ //
// ------------------------------------------------------------------------------ //

// Center and scale the columns to have mean 0 and standard deviation 1
void Normalize(Matrix& mat) {
	int nrows = mat.getNumRows(), ncols = mat.getNumCols();
	
	for (int col = 0; col < ncols; col++) {
		Vector columnVec = mat.getCol(col);
		

		// Compute col mean
		double colMean = mean(columnVec);
		double colStd = stanDev(columnVec, colMean); 

		if (abs(colStd) < 1e-10) { // threshold for saying column is constant
			// only center to 0, don't scale constant column
			for (int coord = 0; coord < nrows; coord++) {
				mat.values[coord][col] = mat.values[coord][col] - colMean; 
			}
		}

		else { // both center and scale column
			for (int coord = 0; coord < nrows; coord++) {
				mat.values[coord][col] = (mat.values[coord][col] - colMean) / colStd;
			}
		}
	}
}

void Normalize(Matrix& mat, Vector& colMeans, Vector& colStds) {
	int nrows = mat.getNumRows(), ncols = mat.getNumCols();

	for (int col = 0; col < ncols; col++) {
		Vector columnVec = mat.getCol(col);

		// Compute col mean
		double colMean = colMeans.values[col];
		double colStd = colStds.values[col];

		if (abs(colStd) < 1e-10) { // threshold for saying column is constant
			for (int coord = 0; coord < nrows; coord++) {
				mat.values[coord][col] = mat.values[coord][col] - colMean;
			}
		}

		else {
			for (int coord = 0; coord < nrows; coord++) {
				mat.values[coord][col] = (mat.values[coord][col] - colMean) / colStd;
			}
		}
	}
}

inline double softThresh(double x, double shrinkage) {
	return sign(x) * max(abs(x) - shrinkage, 0);
}

inline double max(double val1, double val2) {
	if (val1 >= val2)
		return val1;
	else
		return val2;
}