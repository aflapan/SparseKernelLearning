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
		try {
			if(dim < 0)
				throw dim >= 0;
		}
		catch (bool) {
			cout << "Exception -- Vector dimension is negative: " << dim << "\n";
		}

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

	// Swap and assignment operator 
	friend void swap(Vector& vec1, Vector& vec2);

	Vector& operator=(Vector otherVec);

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
	friend double operator*(Vector& vec1, Vector& vec2); // dot product 
	friend Vector operator-(Vector& vec, double scalar); // broadcasting scalar subtraction 
	friend Vector operator+(Vector& vec, double scalar); // broadcasting addition

	// Custom inserter
	friend ostream& operator<<(ostream& stream, Vector& vec);
};


Vector basisVector(int dimension, int coord);

Vector zeroVec(int dimension);