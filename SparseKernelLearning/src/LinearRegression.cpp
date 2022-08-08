/*
A class and functionality for linear models.
The solver methods allow for both Ridge Regression and the Lasso/elastic net.
Uses the Matrix and Vector classes from VectorAndMatrixClasses.cpp
*/

#include <iostream>
#include "Vector.h"
#include "Matrix.h"
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

