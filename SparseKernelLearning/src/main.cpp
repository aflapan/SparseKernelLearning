
#include "Vector.h"
#include "Matrix.h"
#include "MatrixDecompositions.h"
#include "LinearRegression.h"
#include "Kernel.h"
#include <iostream>
using namespace std;


int main()
{
	Matrix trainSamples = readTxt("C:/Users/Alexander/OneDrive//Documents/ProgrammingPractice/Cpp/LinearAlgebra/sample_data.txt");
	int p = trainSamples.getNumCols();
	
	
	Matrix matResponse = readTxt("C:/Users/Alexander/OneDrive/Documents/ProgrammingPractice/Cpp/LinearAlgebra/sample_response.txt");
	Vector response = matResponse.getCol(0); // convert from matrix to vector


	GaussianKernel gk(p, 10.0);

	// Test kernel regression 
	cout << "Forming kernel matrix... ";
	Matrix kernelMat = gk.eval(trainSamples);
	cout << " Done!\n";

	cout << "Centering rows and columsn... ";
	ColRowMeans crMeans = centerRowsCols(kernelMat);
	cout << "Done!\n";

	return 0;
}