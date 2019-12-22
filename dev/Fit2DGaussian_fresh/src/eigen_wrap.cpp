//Eigen wrappers of existing c code for fit2DGaussian function
//Eigen code will facilitate wrapping by pybind to numpy arrays
//Author: Nicolaas van der Voort, AG Seidel, HHU Düsseldorf
//Created December 2, 2019

#include "eigen_wrap.h"
#include "fit2DGaussian.h"
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

/*Eigen::ArrayXd Fit2DGauss_Eigen(
	Eigen::ArrayXd params,
	Eigen::MatrixXd image){
	LVDoubleArray * subimage, M;
	int length;
	MGParam * p;
	double * vars;
	//build needed MGParam struct
	length = image.size();
	subimage = {length, *image.data() };
	M = {length, *image.data() }; //fill with dummy values
	p = { subimage, length, M };
	double * vars = params.data();
	fit2DGaussian(vars, p);
	return Eigen::Map<ArrayXd>(vars);//map vars to Eigen vector
}*/

/*
Eigen::VectorXd fitNew(const Eigen::MatrixXi &img)
{
	Eigen::VectroXd fit;
	
	return fit;
}
*/