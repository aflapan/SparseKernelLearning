#include "Vector.h"
#include "Matrix.h"
#include "LinearRegression.h"
#include "MatrixDecompositions.h"
#include <iostream>


using namespace std;


int main()
{
	Matrix trainSamples = readTxt("C:/Users/Alexander/OneDrive//Documents/ProgrammingPractice/Cpp/LinearAlgebra/sample_data.txt");
	
	Matrix matResponse = readTxt("C:/Users/Alexander/OneDrive/Documents/ProgrammingPractice/Cpp/LinearAlgebra/sample_response.txt");

	Vector response = matResponse.getCol(0); // convert from matrix to vector

	LinearModel model(trainSamples, response);

	LinearSolutions solution = model.train(0.01, 0.2);
	cout << solution;

	return 0;
}