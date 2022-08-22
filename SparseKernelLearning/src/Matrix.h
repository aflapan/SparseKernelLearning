#pragma once

#include <iostream>
#include "Vector.h"
#include <vector>

using namespace std;

// -------------------------------------------------------------------------- //
// ------------------------------ Matrix Class ------------------------------ //
// -------------------------------------------------------------------------- //

class Matrix {
	int nrows;
	int ncols;
public:
	double** values;

	// Constructor and Destructor 
	Matrix() { nrows = ncols = 0; } // Default constructor

	Matrix(int r, int c) {
		nrows = r; // Put exception-handling later
		ncols = c;

		values = new double* [nrows];

		register int row;
		for (row = 0; row < nrows; row++)
			values[row] = new double[ncols];
	}

	~Matrix() {

		register int row;
		for (row = 0; row < nrows; row++) {
			delete[] values[row];
		}
		delete[] values;
	}

	// Copy Constructor 
	Matrix(const Matrix& mat) {
		nrows = mat.nrows, ncols = mat.ncols;

		// initialize memory 
		values = new double* [nrows];

		register int row, col;
		for (row = 0; row < nrows; row++) {
			values[row] = new double[ncols];
		}

		// copy values 
		for (row = 0; row < nrows; row++) {
			for (col = 0; col < ncols; col++) {
				values[row][col] = mat.values[row][col];
			}
		}

	}

	// Swap and assignment operator 
	//friend void swap(Matrix& mat1, Matrix& mat2);

	Matrix& operator=(Matrix otherMat);

	// Functionality to replace a row or column values
	void replaceRow(int rowIndex, Vector& replacement);
	void replaceCol(int colIndex, Vector& replacement);

	// Define summation and scalar-multiplicative operators 
	friend Matrix operator+(Matrix& mat1, Matrix& mat2);
	friend Matrix operator-(Matrix& mat1, Matrix& mat2);
	friend Matrix operator*(Matrix& mat, double scalar);
	friend Matrix operator*(double scalar, Matrix& mat);
	friend Matrix operator/(Matrix& mat, double scalar);

	// transpose function
	friend Matrix transpose(Matrix& mat);

	// matrix multiplication
	friend Vector operator*(Matrix& mat, Vector& vec);
	friend Matrix operator*(Matrix& mat1, Matrix& mat2);

	// Getter functions
	int getNumRows() { return nrows; }
	int getNumCols() { return ncols; }

	// Get individual rows and columns 
	Vector getRow(int rowIndex);
	Vector getCol(int colIndex);

	// Sub-matrix
	Matrix subset(int* rows, int* cols, int numRows, int numCols);

	// Norm functions 
	friend double norm(Matrix& mat);

	// Column Means and Standard Deviations 
	Vector colMeans();
	Vector rowMeans();

	Vector colStanDevs();
	Vector colStanDevs(Vector& colMeans);

	// Custom inserter
	friend ostream& operator<<(ostream& stream, Matrix& mat);

};


// --------------------------------------------------------------------------------------------------------- //
// ------------------------------ Data class for storing row and column means ------------------------------ //
// --------------------------------------------------------------------------------------------------------- //

class ColRowMeans {
	Vector colMeans;
	Vector rowMeans;
public:
	ColRowMeans(Vector cMeans, Vector rMeans) {
		colMeans = cMeans;
		rowMeans = rMeans;
	}

	// getter and setter functions 
	Vector getRowMeans() { return rowMeans; }
	Vector getColMeans() { return colMeans; }

	void setRowMeans(Vector newRowMeans) { rowMeans = newRowMeans; }
	void setColMeans(Vector newColMeans) { colMeans = newColMeans; }
};


// ------------------------------------------------------------------------------------- //
// ------------------------------ Functions and Operators ------------------------------ //
// ------------------------------------------------------------------------------------- //


// Function for centering both the rows and columns. Alters original matrix.
// Returns a data class containing original row and column means. 
ColRowMeans centerRowsCols(Matrix& mat);

// Reading File Input
Matrix readTxt(char const* filename);

// Function to generate the identity matrix 
Matrix Identity(int dimension);

Matrix zeroMat(int nrows, int ncols);

// Fast computation of XtX
Matrix crossprod(Matrix& mat);