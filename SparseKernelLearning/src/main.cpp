
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


	GaussianKernel gk(10.0);
	Vector testVec = trainSamples.getRow(0);
	//Vector kernVec = gk.eval(trainSamples, testVec);

	cout << "Forming kernel matrix... ";
	Matrix kernMat = gk.eval(trainSamples);
	cout << "Done forming kernel matrix!" << "\n";

	return 0;
}