
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


	GaussianKernel gk(p, 100.0);
	KernelRegression kernReg(trainSamples, response, &gk);
	double ridgePen = 1e-10; 
	double sparsityPen = 0.1;
	cout << "Training kernel regression... ";
	kernReg.train(ridgePen, sparsityPen);
	cout << "Done!\n";

	Vector weights = kernReg.getSolution().getWeights();
	cout << weights;
	
}