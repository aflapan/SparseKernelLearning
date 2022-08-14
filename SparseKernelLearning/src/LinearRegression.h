#pragma once


#include <iostream>
#include "Vector.h"
#include "Matrix.h"
using namespace std;

// ------------------------------------------------------------------------------------ //
// ------------------------------ Linear Solutions Class ------------------------------ //
// ------------------------------------------------------------------------------------ //

class LinearSolutions {
	double Intercept;
	Vector Coefficients;
public:
	// Constructor 
	LinearSolutions(double newIntercept, Vector regCoeffs) {
		Intercept = newIntercept;
		Coefficients = regCoeffs;
	}

	 // Getter and setter functions 
	double getIntercept() { return Intercept; }
	Vector getCoefficients() { return Coefficients; }

	void setIntercept(double newIntercept) { Intercept = newIntercept; }
	void setCoefficients(Vector& newCoefficients) { Coefficients = newCoefficients; }
	void updateCoeffVal(double newVal, int coordinate) { Coefficients.values[coordinate] = newVal; }

	// Custom inserter
	friend ostream& operator<<(ostream& stream, LinearSolutions& solution);
};


// -------------------------------------------------------------------------------- //
// ------------------------------ Linear Model Class ------------------------------ //
// -------------------------------------------------------------------------------- //

class LinearModel {
	Matrix trainSamples;
	Vector trainResponse;

public:
	LinearModel(Matrix trainingSamples, Vector trainingResponse) {
		trainSamples = trainingSamples;
		trainResponse = trainingResponse;
	}

	LinearSolutions train();
	LinearSolutions train(double ridgePen);
	LinearSolutions train(double ridgePen, double sparsityPen);

	double objectiveValue(LinearSolutions& intAndCoeffs, double ridgePen, double sparsityPen) {
		int nrows = trainSamples.getNumRows(), ncols = trainSamples.getNumCols();

		// initialize objective function components 
		double ridgeTerm = 0;
		double sparsityTerm = 0;

		// Penalty terms 
		for (int coeff = 0; coeff < ncols; coeff++) {
			double coeffVal = intAndCoeffs.getCoefficients().values[coeff];
			ridgeTerm += coeffVal * coeffVal;
			sparsityTerm += abs(coeffVal);
		}

		// Sum-of-sqaured error
		Vector coefficients = intAndCoeffs.getCoefficients();
		double intercept = intAndCoeffs.getIntercept();

		Vector Xbeta = trainSamples * coefficients;
		Vector linearPrediction = Xbeta + intercept;
		Vector resid = trainResponse - linearPrediction;
		double sumSquaredError = resid * resid;

		return sumSquaredError / nrows + ridgePen * ridgeTerm + sparsityPen * sparsityTerm;
	}
};




// -------------------------------------------------------------------------------- //
// ------------------------------ Solution Functions ------------------------------ //
// -------------------------------------------------------------------------------- //

Vector backSolve(Matrix& upperTriangMat, Vector& targetVec);




// ------------------------------------------------------------------------------ //
// ------------------------------ Helper Functions ------------------------------ //
// ------------------------------------------------------------------------------ //

// Center and scale the columns to have mean 0 and standard deviation 1
void Normalize(Matrix& mat);
void Normalize(Matrix& mat, Vector& colMeans, Vector& colStds);

inline double softThresh(double x, double shrinkage);

inline double max(double val1, double val2);