#pragma once


#include <iostream>
#include "Vector.h"
#include "Matrix.h"
using namespace std;

// -------------------------------------------------------------------------------- //
// ------------------------------ Linear Model Class ------------------------------ //
// -------------------------------------------------------------------------------- //

class LinearModel {
	Matrix TrainingSamples;
	Vector TrainingResponse;

public:
	LinearModel(Matrix trainingSamples, Vector trainingResponse) {
		TrainingSamples = trainingSamples;
		TrainingResponse = trainingResponse;
	}

	Vector train();
	Vector train(double ridgePen);
	Vector train(double ridgePen, double sparsityPen);

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

