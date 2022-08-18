#include <iostream>
#include <cmath>
#include "Vector.h"
#include "Matrix.h"
#include "Kernel.h"
using namespace std;




// ---------------------------------------------------------------------------------------- //
// ------------------------------ Kernel Functional Hash-map ------------------------------ //
// ---------------------------------------------------------------------------------------- //


double gaussianKernel(Vector& data1, Vector& data2, double sigma) {
	Vector diff = data1 - data2;
	double normVal = norm(diff);
	double sumOfSquares = normVal * normVal;
	return exp(-sumOfSquares / sigma)
}

double polynomialKernel(Vector& data1, Vector& data2, int degree) {
	double dotProd = data1 * data2;
	return pow(1 + dotProd, degree);
}





