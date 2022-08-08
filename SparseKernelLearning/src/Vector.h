#pragma once

#include <iostream>

using namespace std;

// -------------------------------------------------------------------------- //
// ------------------------------ Vector Class ------------------------------ //
// -------------------------------------------------------------------------- //


class Vector {
	int dimension;
public:
	double* values;

	// Constructor and Destructor 
	Vector() {dimension = 0;} // default constructor

	Vector(int dim) {
		dimension = dim;
		values = new double[dimension];
	}

	~Vector() {
		delete[] values;
	}

	// getter function for dimension
	int getDim() { return dimension; }

	// Copy Constructor 
	Vector(const Vector& vec) {
		dimension = vec.dimension;
		values = new double[dimension];

		register int coord;
		for (coord = 0; coord < dimension; coord++) {
			values[coord] = vec.values[coord];
		}
	}

	// Euclidean Norm
	double friend norm(Vector& vec);

	// Mean and standard deviation of vectors 
	double friend mean(Vector& vec);

	double friend stanDev(Vector& vec);
	double friend stanDev(Vector& vec, double mean);

	// Override operators 
	friend Vector operator*(Vector& vec, double scalar);
	friend Vector operator*(double scalar, Vector& vec);
	friend Vector operator/(Vector& vec, double scalar);
	friend Vector operator+(Vector& vec1, Vector& vec2);
	friend Vector operator-(Vector& vec1, Vector& vec2);
	friend double operator*(Vector& vec1, Vector& vec2);

	// Custom inserter
	friend ostream& operator<<(ostream& stream, Vector& vec);
};


