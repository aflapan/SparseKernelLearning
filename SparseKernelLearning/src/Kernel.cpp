#include <iostream>
#include <cmath>
#include <stdio.h>

#include "Vector.h"
#include "Matrix.h"
#include "Kernel.h"

using namespace std;




// ------------------------------------------------------------------------------ //
// ------------------------------ Specific Kernels ------------------------------ //
// ------------------------------------------------------------------------------ //



// Gaussian kernel derivative functions 
Vector GaussianKernel::deriv(Vector& vec1, Vector& vec2) {
	Vector diff = vec1 - vec2;
	Vector weights = this->getWeights();
	int dim = this->getDim();

	// Weight difference vectors 
	for (int coord = 0; coord < dim; coord++) {
		diff.values[coord] *= weights.values[coord];
	}

	double multFactor = -(2.0 / scale) * this->eval(vec1, vec2);
	return multFactor * diff;
}

Matrix GaussianKernel::deriv(Matrix& mat, Vector& vec) {
	int nrows = mat.getNumRows(), ncols = mat.getNumCols();
	Matrix derivMat(nrows, ncols);

	Vector weights = this->getWeights();
	int dim = this->getDim();

	for (int row = 0; row < nrows; row++) {
		Vector rowVec = mat.getRow(row);
		Vector diff = rowVec - vec;

		for (int coord = 0; coord < dim; coord++) {
			diff.values[coord] *= weights.values[coord];
		}

		double multFactor = -(2.0 / scale) * this->eval(rowVec, vec);
		diff = multFactor * diff;
		derivMat.replaceRow(row, diff);
	}


	return derivMat;
}



// Gaussian Kernel evaluation functions 

double GaussianKernel::eval(Vector& vec1, Vector& vec2) {
	Vector diff = vec1 - vec2;
	// Scale coefficients of diff vector

	for (int coord = 0; coord < diff.getDim(); coord++) {
		diff.values[coord] *= this->getWeights().values[coord];
	}

	double normVal = norm(diff);
	return exp(-normVal * normVal / scale);
}

Vector GaussianKernel::eval(Matrix& mat, Vector& vec) {
	int dim = vec.getDim();
	int nrows = mat.getNumRows();
	Vector weights = this->getWeights();

	if (mat.getNumCols() != dim)
		throw invalid_argument("Number of matrix columns and vector dimension must be equal.");

	Vector kernelVec(nrows);

	double sqNormVal, diff;

	for (int row = 0; row < nrows; row++) {
		sqNormVal = 0;

		for (int coord = 0; coord < dim; coord++) {
			diff = (vec.values[coord] - mat.values[row][coord]) * weights.values[coord];
			sqNormVal += diff * diff;
		}
		kernelVec.values[row] = exp(-sqNormVal / scale);
	}

	return kernelVec;
}

Matrix GaussianKernel::eval(Matrix& mat1, Matrix& mat2) {
	int ncols1 = mat1.getNumCols(), nrows1 = mat1.getNumRows();
	int ncols2 = mat2.getNumCols(), nrows2 = mat2.getNumRows();
	Vector weights = this->getWeights();

	if (ncols1 != ncols2)
		throw invalid_argument("Number of columns in both matrices must be equal.");

	Matrix kernMat(nrows1, nrows2);

	double sqNormVal, diff;
	for (int row1 = 0; row1 < nrows1; row1++) {
		for (int row2 = 0; row2 < nrows2; row2++) {
			sqNormVal = 0;

			for (int coord = 0; coord < ncols1; coord++) {
				diff = (mat1.values[row1][coord] - mat2.values[row2][coord]) * weights.values[coord];
				sqNormVal += diff * diff;
			}
			kernMat.values[row1][row2] = exp(-sqNormVal / scale);
		}

	}
	return kernMat;
	 
}

Matrix GaussianKernel::eval(Matrix& mat) {
	return eval(mat, mat);
}


// Linear Kernel 

// evaluation functions 
double LinearKernel::eval(Vector& vec1, Vector& vec2) {

	// test dimension correctness 
	int dim1 = vec1.getDim(), dim2 = vec2.getDim(), kernDim = this->getDim();
	if (dim1 != kernDim || dim2 != kernDim)
		throw invalid_argument("Dimension of vectors must be equal to weighted kernel dimension.");

	Vector weights = this->getWeights();
	double weightedDotProd = 0;

	for (int coord = 0; coord < kernDim; coord++) {
		weightedDotProd += weights.values[coord] * weights.values[coord] * vec1.values[coord] * vec2.values[coord];
	}
	return weightedDotProd;
}

Vector LinearKernel::eval(Matrix& mat, Vector& vec) {
	// test dimension correctness 
	int matDim = mat.getNumCols(), vecDim = vec.getDim(), kernDim = this->getDim();
	if (matDim != kernDim || vecDim != kernDim)
		throw invalid_argument("Matrix and Vector must have same dimensionality as weighted kernel.");

	int nrows = mat.getNumRows();
	Vector kernelVec(mat.getNumRows());
	Vector weights = this->getWeights(); 

	for (int row = 0; row < nrows; row++) {
		kernelVec.values[row] = 0;
		for (int col = 0; col < kernDim; col++) {
			kernelVec.values[row] += weights.values[col] * weights.values[col] * mat.values[row][col] * vec.values[col];
		}
	}
	return kernelVec;
}

Matrix LinearKernel::eval(Matrix& mat1, Matrix& mat2) {
	// test dimension correctness 
	int matDim1 = mat1.getNumCols(), mat2Dim = mat2.getNumCols(), kernDim = this->getDim();
	if (matDim1 != kernDim || mat2Dim != kernDim)
		throw invalid_argument("Matrices must have same dimensionality as weighted kernel.");

	int nrows1 = mat1.getNumRows(), nrows2 = mat2.getNumRows();
	Matrix kernelMatrix = zeroMat(nrows1, nrows2);
	Vector weights = this->getWeights();

	for (int row1 = 0; row1 < nrows1; row1++) {
		for (int row2 = 0; row2 < nrows2; row2++) {
			for (int coord = 0; coord < kernDim; coord++) {
				kernelMatrix.values[row1][row2] += weights.values[coord] * weights.values[coord] * mat1.values[row1][coord] * mat2.values[row2][coord];
			}
		}
	}

	return kernelMatrix;
}

Matrix LinearKernel::eval(Matrix& mat) {
	return this->eval(mat, mat);
}


// derivative functions
Vector LinearKernel::deriv(Vector& vec1, Vector& vec2) {

	// test dimension correctness 
	int dim1 = vec1.getDim(), dim2 = vec2.getDim(), kernDim = this->getDim();
	if (dim1 != kernDim || dim2 != kernDim)
		throw invalid_argument("Dimension of vectors must be equal to weighted kernel dimension.");

	Vector kernDeriv(kernDim);
	Vector weights = this->getWeights();

	for (int coord = 0; coord < kernDim; coord++) {
		kernDeriv.values[coord] = 2 * weights.values[coord] * vec1.values[coord] * vec2.values[coord];
	}
	return kernDeriv;
}


Matrix LinearKernel::deriv(Matrix& mat, Vector& vec) {
	// test dimension correctness 
	int matDim = mat.getNumCols(), vecDim = vec.getDim(), kernDim = this->getDim();
	if (matDim != kernDim || vecDim != kernDim)
		throw invalid_argument("Matrix and Vector must have same dimensionality as weighted kernel.");

	int nrows = mat.getNumRows();
	Matrix kernDeriv(nrows, kernDim);
	Vector weights = this->getWeights();

	for (int row = 0; row < nrows; row++) {
		for (int coord = 0; coord < kernDim; coord++) {
			kernDeriv.values[row][coord] = 2 * weights.values[coord] * mat.values[row][coord] * vec.values[coord];
		}
	}
	return kernDeriv;
}

