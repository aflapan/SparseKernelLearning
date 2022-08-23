
#include "Vector.h"
#include "Matrix.h"
#include "MatrixDecompositions.h"
#include "LinearRegression.h"
#include "Kernel.h"
#include "KernelRegression.h"
#include <iostream>
using namespace std;


int main()
{
	Matrix trainSamples = readTxt("C:/Users/Alexander/OneDrive//Documents/ProgrammingPractice/Cpp/LinearAlgebra/sample_data.txt");
	int p = trainSamples.getNumCols(), numSamps = trainSamples.getNumRows();
	
	
	Matrix matResponse = readTxt("C:/Users/Alexander/OneDrive/Documents/ProgrammingPractice/Cpp/LinearAlgebra/sample_response.txt");
	Vector response = matResponse.getCol(0); // convert from matrix to vector


	// Center and scale the training samples and response
	Normalize(trainSamples);
	double responseMean = mean(response);
	response = response - responseMean; 



	GaussianKernel gk(p, 100.0);
	KernelRegression kernReg(trainSamples, response, &gk);
	double ridgePen = 1e-10; 
	double sparsityPen = 1e-4;

	kernReg.train(ridgePen, sparsityPen);
	Vector finalWeights = kernReg.getSolution().getWeights();

	cout << finalWeights << "\n\n";
	
	/*
	Matrix testQ = Identity(p);
	Vector testBeta = zeroVec(p);
	Vector initialWeights = zeroVec(p);

	for (int coord = 0; coord < p; coord++) {
		double val = 1 / ((double) coord + 1);
		testQ.values[coord][coord] = val;
		testBeta.values[coord] = 2*val;
	}
	
	cout << "\nMinimizing quadratic form... ";

	Vector testWeights = minimizeQuadForm(testQ, testBeta, initialWeights, 0.0);
	cout << testWeights;
	*/

}