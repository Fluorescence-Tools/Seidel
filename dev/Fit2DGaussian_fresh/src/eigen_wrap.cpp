//Eigen wrappers of existing c code for fit2DGaussian function
//Eigen code will facilitate wrapping by pybind to numpy arrays
//Author: Nicolaas van der Voort, AG Seidel, HHU Düsseldorf
//Created December 2, 2019

#include "eigen_wrap.h"
#include <Eigen/Core>
#include <Eigen/LU>


using namespace Eigen;

// takes numpy array as input and returns another numpy array
Eigen::MatrixXd inv(Eigen::MatrixXd xs) {
	return xs.inverse();
}

int subtract(int i, int j)
{
	return i - j;
}

/*
Eigen::VectorXd fitNew(const Eigen::MatrixXi &img)
{
	Eigen::VectroXd fit;
	
	return fit;
}
*/