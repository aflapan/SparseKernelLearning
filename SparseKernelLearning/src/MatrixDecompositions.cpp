/*
This script contains functionality for the common matrix decompositions,
such as reduced QR, LU, SVD, etc.

The implementations are based on the reference 
"Numerical Linear Algebra", Lloyd N. Trefethen and David Bau III, SIAM, 1997.
*/

#include "Matrix.h"
#include "Vector.h"
#include "MatrixDecompositions.h"

#include <iostream>
using namespace std;



// ------------------------------------------------------------------------------ //
// ------------------------------ QR Decomposition ------------------------------ //
// ------------------------------------------------------------------------------ //



/*
Implements the Modified Gram-Schmidt Algorithm for reduced-form QR decomposition.

Note, overwrites original matrix to have orthonormal column vectors, transforming
it into the matrix Q. Returns the upper-triangular square matrix R satisfying QR = original matrix. 
*/

QR qrDecomp(Matrix mat) {
	int nrows = mat.getNumRows(), ncols = mat.getNumCols();
	Matrix R(ncols, ncols);

	register int col, nextCol, row;

	for (col = 0; col < ncols; col++) {
		double sqColNorm = 0;

		// Compute norms 
		for (row = 0; row < nrows; row++) {
			sqColNorm += mat.values[row][col] * mat.values[row][col];
		}

		R.values[col][col] = pow(sqColNorm, 0.5);

		// Scale column by norm value
		for (row = 0; row < nrows; row++) {
			mat.values[row][col] = mat.values[row][col] / R.values[col][col];
		}


		// subtract projections onto Q column from subsequent columns in mat.
		for (nextCol = col + 1; nextCol < ncols; nextCol++) {

			double dotprod = 0;
			register int rowIndex;

			for (rowIndex = 0; rowIndex < nrows; rowIndex++) {
				dotprod = mat.values[rowIndex][col] * mat.values[rowIndex][nextCol];
			}

			R.values[col][nextCol] = dotprod;

			for (rowIndex = 0; rowIndex < nrows; rowIndex++) {
				mat.values[rowIndex][nextCol] = mat.values[rowIndex][nextCol] - R.values[col][nextCol] * mat.values[rowIndex][col];
			}

		}
	}

	QR qr(mat, R);

	return qr;
}

Matrix Householder(Matrix mat) {

	int ncols = mat.getNumCols(), nrows = mat.getNumRows();

	Matrix ReflectionVectors = zeroMat(nrows, ncols);

	for (int col = 0; col < ncols; col++) {
		int subDim = nrows - col;

		Vector colSubVec = zeroVec(nrows);

		for (int row = col; row < nrows; row++) {
			colSubVec.values[row] = mat.values[row][col]; 
		}

		Vector basisVec = basisVector(nrows, col);
		double scalar = sign(colSubVec.values[0]) * norm(colSubVec);
		
		basisVec = scalar * basisVec;
		Vector orthogVec = basisVec + colSubVec;
		orthogVec = orthogVec / norm(orthogVec);
		
		// Save orthonormal columns 
		ReflectionVectors.replaceCol(col, orthogVec);

		for (int nextCol = col; nextCol < ncols; nextCol++) {

			Vector nextColSubVec = zeroVec(nrows);
			for (int row = col; row < nrows; row++) {
				nextColSubVec.values[row] = mat.values[row][nextCol];
			}

			Vector reflection = orthogVec * (2 * (orthogVec * nextColSubVec));
			reflection = nextColSubVec - reflection; // issue here
													 
			// replace values in matrix
			for (int row = col; row < nrows; row++) {
				mat.values[row][nextCol] = reflection.values[row-col]; 
			}
		}
	

	}

	return ReflectionVectors;

}




// ------------------------------------------------------------------------------ //
// ------------------------------ Helper Functions ------------------------------ //
// ------------------------------------------------------------------------------ //



inline int sign(double x) {
	if (x >= 0)
		return 1;
	else
		return -1;
}