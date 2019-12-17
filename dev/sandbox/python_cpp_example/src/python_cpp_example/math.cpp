#include <Eigen/Core>
#include <Eigen/LU>
#include <iostream>
#include "math.hpp"


using namespace std;
using namespace Eigen;

int add(int i, int j)
{
	return i + j;
}

int subtract(int i, int j)
{
	return i - j;
}

// takes numpy array as input and returns another numpy array
Eigen::MatrixXd inv(Eigen::MatrixXd xs) {
	return xs.inverse();
}