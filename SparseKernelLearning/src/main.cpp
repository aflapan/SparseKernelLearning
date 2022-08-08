#include "Vector.h"
#include "Matrix.h"
#include "LinearRegression.h"
#include "MatrixDecompositions.h"
#include <iostream>


using namespace std;


int main()
{
	Matrix trainSamples = readTxt("C:/Users/Alexander/OneDrive//Documents/ProgrammingPractice/Cpp/LinearAlgebra/sample_data.txt");
	Matrix response = readTxt("C:/Users/Alexander/OneDrive/Documents/ProgrammingPractice/Cpp/LinearAlgebra/sample_response.txt");

	Vector vecResponse = response.getCol(0); // convert from matrix to column
	Matrix R = QR(trainSamples);
	Matrix transposeSamples = transpose(trainSamples);
	Vector targetVec = transposeSamples * vecResponse;
	Vector solutions = backSolve(R, targetVec);

	cout << solutions;

	return 0;
}