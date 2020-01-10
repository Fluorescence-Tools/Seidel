//Eigen wrappers of existing c code for fit2DGaussian function
//Eigen code will facilitate wrapping by pybind to numpy arrays
//Author: Nicolaas van der Voort, AG Seidel, HHU Düsseldorf
//Created December 2, 2019

#include "pywrap.h"
#include "fit2DGaussian.h"
#include "twoIstar.h"
#include <Eigen/Core>
#include <Eigen/LU>
#include <iostream>


using namespace Eigen;

// takes numpy array as input and returns another numpy array
Eigen::MatrixXd inv(Eigen::MatrixXd xs) {
	return xs.inverse();
}

int subtract(int i, int j)
{
	return (i - j);
}

/*wrapper for Fit2DGaussian function.
params must contain 18(!) parameters, and can be called from python using a 1D
numpy array.
data contains the image to be fitted and can be called from python using
a 2D numpy array.
No new memory is created in this wrapper. Any memory allocated by subroutines
is freed by subroutines.
returns optimsed params vector*/
ArrayXd Fit2DGaussian_pywrap(
	ArrayXd vars,
	MatrixXd data
) {
	//(ncols, nrows) -> (xlen, ylen)
	//be mindful of a mixup.
	//if correct,	vars[0,6,9] represents x
	//				vars[1,7,10] represents y
	double* vars_p;
	double* data_p;
	int ncols;
	int nrows;

	//get pointer to data
	ncols = (int)data.cols();
	nrows = (int)data.rows();
	vars_p = vars.data();
	data_p = data.data();

	//if mixup, replace with fit2DGaussian(vars_p, data_p, ncols, nrows);
	fit2DGaussian(vars_p, data_p, nrows, ncols);

	return vars;
}

/*
params must contain 6 values to model a single Gaussian
returns a filled model.
*/
MatrixXd model2DGaussian_pywrap(
	ArrayXd vars,
	MatrixXd model
) {
	double* vars_p;
	double* model_p;
	int ncols;
	int nrows;

	//get pointer to data
	ncols = (int)model.cols();
	nrows = (int)model.rows();
	vars_p = vars.data();
	model_p = model.data();

	//call function
	model2DGaussian(vars_p, model_p, nrows, ncols);

	return model;
}

/*
params must contain 9 values to model two Gaussians
returns a filled model.
*/
MatrixXd modelTwo2DGaussian_pywrap(
	ArrayXd vars,
	MatrixXd model
) {
	double* vars_p;
	double* model_p;
	int ncols;
	int nrows;

	//get pointer to data
	ncols = (int)model.cols();
	nrows = (int)model.rows();
	vars_p = vars.data();
	model_p = model.data();

	//call function
	modelTwo2DGaussian(vars_p, model_p, nrows, ncols);

	return model;
}

/*
params must contain 12 values to model three Gaussians
returns a filled model.
*/
MatrixXd modelThree2DGaussian_pywrap(
	ArrayXd vars,
	MatrixXd model
) {
	double* vars_p;
	double* model_p;
	int ncols;
	int nrows;

	//get pointer to data
	ncols = (int)model.cols();
	nrows = (int)model.rows();
	vars_p = vars.data();
	model_p = model.data();

	//call function
	modelThree2DGaussian(vars_p, model_p, nrows, ncols);

	return model;
}

/*
get two I star value for data modelled by model.
Note: therms depending on data only have been dropped in the calculation
therefore values may differ from other optimisers.
*/
double W2DG_pywrap(
	MatrixXd data,
	MatrixXd model
) {
	double* data_p;
	double* model_p;
	int ncols;
	int nrows;

	//get pointer to data
	ncols = (int)model.cols();
	nrows = (int)model.rows();
	data_p = data.data();
	model_p = model.data();

	return W2DG(data_p, model_p, nrows * ncols);
}