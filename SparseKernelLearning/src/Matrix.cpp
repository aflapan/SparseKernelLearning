#include "Matrix.h"
#include "Vector.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>

using namespace std;

// custom inserter
ostream& operator<<(ostream& stream, Matrix& mat) {
	int nrows = mat.getNumRows(), ncols = mat.getNumCols();

	for (int row = 0; row < nrows; row++) {
		for (int col = 0; col < ncols; col++) {
			stream << mat.values[row][col] << ' ';
		}
		stream << "\n";
	}
	return stream;
}

/* Define Matrix Operations */

Matrix operator*(Matrix& mat1, Matrix& mat2) {

	int nrow1 = mat1.getNumRows(), ncol1 = mat1.getNumCols();
	int nrow2 = mat2.getNumRows(), ncol2 = mat2.getNumCols();

	Matrix prod(nrow1, ncol2);

	for (int row = 0; row < nrow1; row++) {
		for (int col = 0; col < ncol2; col++) {

			double dotprod = 0;
			for (int intermediate = 0; intermediate < ncol1; intermediate++) {
				dotprod += mat1.values[row][intermediate] * mat2.values[intermediate][col];
			}
			prod.values[row][col] = dotprod;
		}
	}
	return prod;
}

Vector operator*(Matrix& mat, Vector& vec) {
	int ncols = mat.getNumCols(), nrows = mat.getNumRows();
	Vector imageVec(nrows);

	for (int row = 0; row < nrows; row++) {
		double dotProd = 0;
		for (int col = 0; col < ncols; col++) {
			dotProd += mat.values[row][col] * vec.values[col];
		}
		imageVec.values[row] = dotProd;
	}
	return imageVec;
}

Matrix operator+(Matrix& mat1, Matrix& mat2) {
	int nrows = mat1.getNumRows(), ncols = mat1.getNumCols();
	Matrix sum(nrows, ncols);
	// element-wise sum of terms
	register int row, col;
	for (row = 0; row < nrows; row++) {
		for (col = 0; col < ncols; col++) {
			sum.values[row][col] = mat1.values[row][col] + mat2.values[row][col];
		}
	}
	return sum;
}

Matrix operator-(Matrix& mat1, Matrix& mat2) {
	int nrows = mat1.getNumRows(), ncols = mat1.getNumCols();

	Matrix matDiff(nrows, ncols);
	// element-wise sum of terms
	register int row, col;
	for (row = 0; row < nrows; row++) {
		for (col = 0; col < ncols; col++) {
			matDiff.values[row][col] = mat1.values[row][col] - mat2.values[row][col];
		}
	}
	return matDiff;
}

Matrix operator*(Matrix& mat, double scalar) {
	int nrows = mat.getNumRows(), ncols = mat.getNumCols();
	Matrix prod(nrows, ncols);

	register int row, col;
	for (row = 0; row < nrows; row++) {
		for (col = 0; col < ncols; col++) {
			prod.values[row][col] = scalar * mat.values[row][col];
		}
	}

	return prod;
}

Matrix operator*(double scalar, Matrix& mat) {
	int nrows = mat.getNumRows(), ncols = mat.getNumCols();
	Matrix prod(nrows, ncols);

	register int row, col;
	for (row = 0; row < nrows; row++) {
		for (col = 0; col < ncols; col++) {
			prod.values[row][col] = scalar * mat.values[row][col];
		}
	}

	return prod;
}

Matrix operator/(Matrix& mat, double scalar) {
	int nrows = mat.getNumRows(), ncols = mat.getNumCols();
	Matrix dividedMat(nrows, ncols);

	register int row, col;
	for (row = 0; row < nrows; row++) {
		for (col = 0; col < ncols; col++) {
			dividedMat.values[row][col] = mat.values[row][col] / scalar;
		}
	}
	return dividedMat;
}


// Swap and assignment operator 
void swap(Matrix& mat1, Matrix& mat2) {
	using std::swap;

	swap(mat1.nrows, mat2.nrows);
	swap(mat2.ncols, mat2.ncols);
	swap(mat1.values, mat2.values);
}

Matrix& Matrix::operator=(Matrix otherMat) {

	// Delete old values 
	for (int row = 0; row < nrows; row++) {
		delete[] values[row];
	}
	delete[] values;

	// allocate new memory
	nrows = otherMat.getNumRows();
	ncols = otherMat.getNumCols();

	double** values;
	values = new double* [nrows];

	register int row;
	for (row = 0; row < nrows; row++)
		values[row] = new double[ncols];

	// Copy values 
	for (int row = 0; row < nrows; row++) {
		for (int col = 0; col < ncols; col++) {
			values[row][col] = otherMat.values[row][col];
		}
	}

	this->values = values;
	return *this;
}


Matrix transpose(Matrix& mat) {
	int nrows = mat.getNumRows(), ncols = mat.getNumCols();
	Matrix trans(ncols, nrows);

	register int row, col;
	for (col = 0; col < ncols; col++) {
		for (row = 0; row < nrows; row++) {
			trans.values[col][row] = mat.values[row][col];
		}
	}
	return trans;
}

// Fucntionality for replacing a row or column in a matrix 
void Matrix::replaceRow(int rowIndex, Vector& replacement) {
	register int col;

	for (col = 0; col < ncols; col++) {
		values[rowIndex][col] = replacement.values[col];
	}
}

void Matrix::replaceCol(int colIndex, Vector& replacement) {
	register int row;

	for (row = 0; row < nrows; row++) {
		values[row][colIndex] = replacement.values[row];
	}
}


// Subsettting rows and columsn 
Vector Matrix::getRow(int rowIndex) {
	Vector rowVec(ncols);
	register int col;

	for (col = 0; col < ncols; col++) {
		rowVec.values[col] = values[rowIndex][col];
	}
	return rowVec;
}

Vector Matrix::getCol(int colIndex) {
	Vector colVec(nrows);
	register int row;

	for (row = 0; row < nrows; row++) {
		colVec.values[row] = values[row][colIndex];
	}
	return colVec;
}

Matrix Matrix::subset(int* rows, int* cols, int numRows, int numCols) {
	Matrix subMat(numRows, numCols);

	register int row, col;
	for (row = 0; row < numRows; row++) {
		for (col = 0; col < numCols; col++) {
			int rowIndex = rows[row], colIndex = cols[col];
			subMat.values[row][col] = values[rowIndex][colIndex];
		}
	}
	return subMat;
}


// Column Means and Standard Deviations 
Vector Matrix::colMeans() {
	Vector cMeans(ncols);

	for (int colIndex = 0; colIndex < ncols; colIndex++) {
		Vector colVec = this->getCol(colIndex);
		cMeans.values[colIndex] = mean(colVec);
	}
	return cMeans;
}

Vector Matrix::rowMeans() {
	Vector rMeans(nrows);
	for (int rowIndex = 0; rowIndex < nrows; rowIndex++) {
		Vector rowVec = this->getRow(rowIndex);
		rMeans.values[rowIndex] = mean(rowVec);
	}
	return rMeans;
}


Vector Matrix::colStanDevs() {
	Vector stanDevs(ncols);

	for (int colIndex = 0; colIndex < ncols; colIndex++) {
		Vector colVec = this->getCol(colIndex);
		double colMean = mean(colVec);
		stanDevs.values[colIndex] = stanDev(colVec, colMean);
	}
	return stanDevs;
}

Vector Matrix::colStanDevs(Vector& colMeans) {
	Vector stanDevs(ncols);
	for (int colIndex = 0; colIndex < ncols; colIndex++) {
		Vector colVec = this->getCol(colIndex);
		stanDevs.values[colIndex] = stanDev(colVec, colMeans.values[colIndex]);
	}
	return stanDevs;
}


// Frobenius Norm Function 
double norm(Matrix& mat) {
	double squaredNorm = 0.0;
	int nrows = mat.getNumRows(), ncols = mat.getNumCols();

	register int row, col;
	for (row = 0; row < nrows; row++) {
		for (col = 0; col < ncols; col++) {
			squaredNorm += mat.values[row][col] * mat.values[row][col];
		}
	}
	return pow(squaredNorm, 0.5);
}


// Reading File Input
Matrix readTxt(char const* filename) {
	cout << "\nReading data... ";
	int rowCount = 0, colCount = 0, totalSize = 0;

	ifstream fin(filename, ios::in | ios::binary);

	string line;

	int nrows = 0;
	int ncols;

	std::vector<double> valueHolder;
	double val;

	while (getline(fin, line)) {


		std::stringstream ss(line);
		while (ss >> val)
			valueHolder.push_back(val);
		nrows++;
	}

	fin.close();
	ncols = valueHolder.size() / nrows;

	Matrix mat(nrows, ncols);

	for (int row = 0; row < nrows; row++) {
		for (int col = 0; col < ncols; col++) {
			mat.values[row][col] = valueHolder[row * ncols + col];
		}
	}
	cout << "Done!\n";
	return mat;
}

// ----------------------------------------------------------------------------------------- //
// ------------------------------ Additional Matrix Functions ------------------------------ //
// ----------------------------------------------------------------------------------------- //


// Function for centering both the rows and columns. Alters original matrix.
// Returns a data class containing original row and column means. 
ColRowMeans centerRowsCols(Matrix& mat) {
	int nrows = mat.getNumRows(), ncols = mat.getNumCols();

	// Center columns
	Vector cMeans = mat.colMeans();
	for (int colIndex = 0; colIndex < ncols; colIndex++) {
		Vector column = mat.getCol(colIndex);
		column = column - cMeans.values[colIndex];
		mat.replaceCol(colIndex, column);
	}

	// Center rows 
	Vector rMeans = mat.rowMeans();
	for (int rowIndex = 0; rowIndex < nrows; rowIndex++) {
		Vector row = mat.getRow(rowIndex);
		row = row - rMeans.values[rowIndex];
		mat.replaceRow(rowIndex, row);
	}
	// store original row and column values 
	ColRowMeans colAndRowMeans(cMeans, rMeans);
	return colAndRowMeans;
}

// Function to generate the identity matrix 
Matrix Identity(int dimension) {
	Matrix identityMat(dimension, dimension);

	// Initialize all entries to 0
	register int row, col;
	for (row = 0; row < dimension; row++) { // not if-guarding row != col to avoid O(n^2) checks 
		for (col = 0; col < dimension; col++) {
			identityMat.values[row][col] = 0.0;
		}
	}

	// Set diagonal to 1.
	for (row = 0; row < dimension; row++) {
		identityMat.values[row][row] = 1.0;
	}

	return identityMat;
}

Matrix zeroMat(int nrows, int ncols) {
	Matrix zero(nrows, ncols);

	for (int row = 0; row < nrows; row++) {
		for (int col = 0; col < ncols; col++) {
			zero.values[row][col] = 0;
		}
	}
	return zero;
}