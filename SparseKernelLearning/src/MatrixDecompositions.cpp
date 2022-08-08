/*
This script contains functionality for the common matrix decompositions,
such as reduced QR, LU, SVD, etc.

The implementations are based on the reference 
"Numerical Linear Algebra", Lloyd N. Trefethen and David Bau III, SIAM, 1997.
*/

#include "Matrix.h"
#include "Vector.h"

#include <iostream>
using namespace std;


/*
Implements the Modified Gram-Schmidt Algorithm for reduced-form QR decomposition.

Note, overwrites original matrix to have orthonormal column vectors, transforming
it into the matrix Q. Returns the upper-triangular square matrix R satisfying QR = original matrix. 
*/

Matrix QR(Matrix& mat) {
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

	return R;
}



