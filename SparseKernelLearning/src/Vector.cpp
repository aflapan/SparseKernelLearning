/*
A vector and matrix class along with functions.
Author: Alexander Lapanowski <aflapan@gmail.com>
*/

#include <iostream>
#include <cmath>

#include "Vector.h"

using namespace std;

// Swap function
void swap(Vector& vec1, Vector& vec2) {
	using std::swap;
	swap(vec1.dimension, vec2.dimension);
	swap(vec1.values, vec2.values);
}

Vector& Vector::operator=(Vector otherVec) {
	// deallocate memory
	delete[] values;

	// rellocate memory and values 
	dimension = otherVec.getDim();

	double* values;
	values = new double[dimension];

	for (int coord = 0; coord < dimension; coord++) {
		values[coord] = otherVec.values[coord];
	}

	this->values = values;
	return *this;
}


double norm(Vector& vec) {
	double squaredNorm = 0.0;
	int dim = vec.getDim();

	register int coord;
	for (coord = 0; coord < dim; coord++) {
		squaredNorm += vec.values[coord] * vec.values[coord];
	}
	return pow(squaredNorm, 0.5);
}

// Overriding operators
Vector operator*(Vector& vec, double scalar) {
	int dim = vec.getDim();
	Vector scaledVec(dim);

	register int coord;
	for (coord = 0; coord < dim; coord++) {
		scaledVec.values[coord] = vec.values[coord] * scalar;
	}
	return scaledVec;
}

Vector operator*(double scalar, Vector& vec) {

	int dim = vec.getDim();
	Vector scaledVec(dim);

	register int coord;
	for (coord = 0; coord < dim; coord++) {
		scaledVec.values[coord] = vec.values[coord] * scalar;
	}
	return scaledVec;
}

Vector operator/(Vector& vec, double scalar) {
	int dim = vec.getDim();
	Vector scaledVec(dim);

	register int coord;
	for (coord = 0; coord < dim; coord++) {
		scaledVec.values[coord] = vec.values[coord] / scalar;
	}
	return scaledVec;
}

double operator*(Vector& vec1, Vector& vec2) {
	int dim = vec1.getDim();
	double dotprod = 0.0;

	register int i;
	for (i = 0; i < dim; i++) {
		dotprod += vec1.values[i] * vec2.values[i];
	}
	return dotprod;
}

Vector operator+(Vector& vec1, Vector& vec2) {
	int dim = vec1.getDim();
	Vector sumVec(dim);

	register int coord;
	for (coord = 0; coord < dim; coord++) {
		sumVec.values[coord] = vec1.values[coord] + vec2.values[coord];
	}
	return sumVec;
}

Vector operator-(Vector& vec1, Vector& vec2) {
	int dim = vec1.getDim();
	Vector diffVec(dim);

	register int coord;
	for (coord = 0; coord < dim; coord++) {
		diffVec.values[coord] = vec1.values[coord] - vec2.values[coord];
	}
	return diffVec;
}

Vector operator-(Vector& vec, double scalar) {
	int dim = vec.getDim();
	Vector shiftedVec(dim);

	for (int i = 0; i < dim; i++) {
		shiftedVec.values[i] = vec.values[i] - scalar;
	}
	return shiftedVec;
}

Vector operator+(Vector& vec, double scalar) {
	int dim = vec.getDim();
	Vector shiftedVec(dim);

	for (int i = 0; i < dim; i++) {
		shiftedVec.values[i] = vec.values[i] + scalar;
	}
	return shiftedVec;
}

// custom inserter
ostream& operator<<(ostream& stream, Vector& vec) {
	int dim = vec.getDim();

	for (int coord = 0; coord < dim; coord++) {
		stream << vec.values[coord] << ' ';
	}
	return stream;
}

// Vector mean and standard deviation
double mean(Vector& vec) {
	int dim = vec.getDim();
	double sum = 0;

	for (int coord = 0; coord < dim; coord++) {
		sum += vec.values[coord];
	}
	return sum / (double)dim;
}

double stanDev(Vector& vec) {
	int dim = vec.getDim();
	double vecMean = mean(vec);

	double sumOfSquares = 0;
	for (int coord = 0; coord < dim; coord++) {
		sumOfSquares += (vec.values[coord] - vecMean) * (vec.values[coord] - vecMean);
	}

	return pow(sumOfSquares / dim, 0.5);
}

double stanDev(Vector& vec, double mean) {
	int dim = vec.getDim();

	double sumOfSquares = 0;
	for (int coord = 0; coord < dim; coord++) {
		sumOfSquares += (vec.values[coord] - mean) * (vec.values[coord] - mean);
	}

	return pow(sumOfSquares / dim, 0.5);
}


// --------------------------------------------------------------------------------------- //
// ------------------------------ Vector Creation Functions ------------------------------ //
// --------------------------------------------------------------------------------------- //



Vector basisVector(int dimension, int coord) {
	Vector basisVec(dimension);

	for (int coord = 0; coord < dimension; coord++) {
		basisVec.values[coord] = 0;
	}

	basisVec.values[coord] = 1;
	return basisVec;
}

Vector zeroVec(int dimension) {
	Vector zero(dimension);
	for (int coord = 0; coord < dimension; coord++) {
		zero.values[coord] = 0;
	}
	return zero;
}