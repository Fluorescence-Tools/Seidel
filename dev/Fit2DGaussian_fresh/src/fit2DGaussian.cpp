// Fit 2DGauss
// Version 2019-04-16
//edited extensively by Nicolaas van der Voort
//AG Seidel, Düsseldorf

#include "fit2DGaussian.h"
#include "i_lbfgs.h"
#include "twoIstar.h"
#include <iostream>
#include <fstream>
#include <math.h>
//#include <set>
#include <Eigen/Core>

using namespace std;

const int NVARS = 6;
const int NPEAKS_FACTOR = 10;

typedef Eigen::Matrix<double, 5, 1> Vector5d;
typedef Eigen::Matrix<double, 6, 1> Vector6d;

// background and threshold surfaces
static double* bg_surface = NULL;
static int* threshold_surface = NULL;
static int last_threshold = -108;

int osize;																		//object size
int osize_sq;																	//object area
int* icut;

void save_array(double *a, int l)
{	
	ofstream myfile("C:\\Temp\\Headers\\Gauss2D.dat");
		if (myfile.is_open())
		{
			for (int i = 0; i < l; i++) {
				myfile << a[i] << " ";
			}
			myfile.close();
		}
		else cout << "Unable to open file";
}


// weighted residuals, weights*(model-data)

/*
///////////////////////////////// 2DGaussian ////////////////////////////////////////////////
void modelGaussian(Eigen::MatrixXd model, Vector5d params)
{
	//fill the matrix model with a Gaussian according to params
	//put parameters in more descriptive wordings
	int x, y;
	const int cols = model.cols();
	double* ex = new double[cols];
	double ey;
	double x0 = params(0);
	double y0 = params(1);
	double A = params(2);
	double sigma = params(3);
	double bg = params(4);
	double tx;
	double ty;
	tx = 0.5 / (sigma * sigma);
	ty = 0.5 / (sigma * sigma *  * vars[4]);

	for (x = 0; x < model.rows(); x++)
		ex[x] = exp(-(x - x0)*(x - x0)*tx);

	for (y = 0; y < model.cols(); y++) {
		ey = exp(-(y - y0) * (y - y0) * ty);
		for (x = 0; x < model.rows(); x++) {
			model(x, y) = A * ex[x] * ey + bg;
		}
	}
	delete[] ex;

}*/

///////////////////////////////// target function (to minimize) ///////////////////////////////
/*
Function to minimize by bfgs object.
bfgs constructor needs function with arguments (double *, void*)
This function resets parameters within bounds, calculates model and gets goodness
Maybe these functionalities should be split up further?
*/
double target2DGaussian(double* vars, void* gdata_dummy)
{
	double w;
	//convert void into GaussDataType. 
	GaussDataType* gdata = (GaussDataType*)gdata_dummy; 
	int osize = gdata->xlen * gdata->ylen;
	//	MGParam* p = (MGParam*)pM;
	//	LVDoubleArray *subimage = *(p->subimage), *M = *(p->M);
	vars[0] = varinbounds(vars[0], 0, (double)gdata->xlen);
	vars[1] = varinbounds(vars[1], 0, (double)gdata->ylen);
	vars[6] = varinbounds(vars[6], 0, (double)gdata->xlen);
	vars[7] = varinbounds(vars[7], 0, (double)gdata->ylen);
	vars[9] = varinbounds(vars[9], 0, (double)gdata->xlen);
	vars[10] = varinbounds(vars[10], 0, (double)gdata->ylen);
	vars[5] = varlowerbound(vars[5], 0); //if bg <0, bg = 1
 
	//get model
	if ((int) vars[16] == 0 ) {
		model2DGaussian(vars, gdata->model, gdata->xlen, gdata->ylen);
	}
	else if ((int) vars[16] == 1) {
		modelTwo2DGaussian(vars, gdata->model, gdata->xlen, gdata->ylen);
	}
	else if ((int) vars[16] == 2) {
		modelThree2DGaussian(vars, gdata->model, gdata->xlen, gdata->ylen);
	}

	w = W2DG(gdata->data, gdata->model, osize);
	return w;
}

/*
check if var is within bounds
If out-of-bound, reset parameter to the middle of the bounds
*/
double varinbounds(double var, double min, double max) {
	if (var < min || var> max)
		var = (max - min) / 2;
	return var;
}
/*
check if var is below lower bound
if yes, reset to lower-bound + 1*/
double varlowerbound(double var, double min) {
	if (var < min)
		var = min + 1;
	return var;
}
//////////////////////////////////////////// model2DGaussian ////////////////////////////////////////////
// vars = [x0 y0 A sigma ellipticity bg]
int model2DGaussian(double* vars, double * model, int xlen, int ylen)
{
	//fill the matrix model with a Gaussian according to params
	//put parameters in more descriptive wordings
	int x, y;
	int i;
	double* ex = new double[xlen];
	double ey;
	double x0 = vars[0];
	double y0 = vars[1];
	double A = vars[2];
	double sigma = vars[3];
	double eps = vars[4];
	double bg = vars[5];
	double tx;
	double ty;
	tx = 0.5 / (sigma * sigma);
	ty = 0.5 / (sigma * sigma * eps * eps);

	//f(x,y) cn be written as ex(x)*ey(y)
	//first calc ex(x)
	for (x = 0; x < xlen; x++)
		ex[x] = exp(-(x - x0) * (x - x0) * tx);
	//now calc whole thing
	i = 0;
	for (y = 0; y < ylen; y++) {
		ey = exp(-(y - y0) * (y - y0) * ty);
		for (x = 0; x < xlen; x++) {
			model[i] = A * ex[x] * ey + bg;
			i++;
		}
	}
	delete[] ex;
	return 0;
}


//////////////////////////////////////////// modelTwo2DGaussian ////////////////////////////////////////////
//function uses model2DGaussian function for constructor, see fit2DGaussian for parameter declaration
int modelTwo2DGaussian(double* vars, double * model, int xlen, int ylen) {
	int osize = xlen * ylen;
	double * vars_dummy = new double[6];
	double * model_dummy = new double[osize];
	int i;
	//create the first Gaussian with bg, store in model
	model2DGaussian(vars, model, xlen, ylen);

	//create the second Gaussian without bg, store in M_dummy
	vars_dummy[0] = vars[6];
	vars_dummy[1] = vars[7];
	vars_dummy[2] = vars[8];
	vars_dummy[3] = vars[3];
	vars_dummy[4] = vars[4];
	vars_dummy[5] = 0;
	model2DGaussian(vars_dummy, model_dummy, xlen, ylen);

	//add
	for (i = 0; i < osize; i++) {
		model[i] += model_dummy[i];
	}

	delete[] vars_dummy, model_dummy;
	return 0;
}

//////////////////////////////////////////// modelThree2DGaussian ////////////////////////////////////////////

//function uses model2DGaussian and model2DGaussian function for constructor, see fit2DGaussian for parameter declaration
int modelThree2DGaussian(double* vars, double* model, int xlen, int ylen) {
	int osize = xlen * ylen;
	double * vars_dummy = new double[6];
	double * model_dummy = new double[osize];
	int i;
	//create the first two Gaussian with bg, store in M
	modelTwo2DGaussian(vars, model, xlen, ylen);

	//create the third Gaussian without bg, store in M_dummy
	vars_dummy[0] = vars[9];
	vars_dummy[1] = vars[10];
	vars_dummy[2] = vars[11];
	vars_dummy[3] = vars[3];
	vars_dummy[4] = vars[4];
	vars_dummy[5] = 0;

	model2DGaussian(vars_dummy, model_dummy, xlen, ylen);

	//add
	for (i = 0; i < osize; i++) {
		model[i] += model_dummy[i];
	}

	delete[] vars_dummy, model_dummy;
	return 0;
}

//////////////////////////////////////////// fit2DGaussian ////////////////////////////////////////////
/*
//fit2DGaussian initializes optimisation routine
//vars contains the parameters that are optimized and has length 18!
//  vars:
	0: x0, 
	1: y0, 
	2: A0, 
	3: sigma, 
	4: ellipticity, 
	5: bg, 
	6: x1, 
	7: y1, 
	8: A1, 
	9: x2, 
	10: y2, 
	11: A2, 
	12: info, contains information from the fitting algorithm
	13: wi_nowi, outdated
	14: fit_bg, asks if background is fitted. 0 -> bg is fitted
	15: ellipt_circ, determines if elliptical fits are allowed. Only relevant for model2DGaussian
	16: model, determines the model to be used:
	    0: model2DGaussian
		1: modelTwo2DGaussian
		2: modelThree2DGaussian
	17: reserved for two Istar value of optimised solution
*/
int fit2DGaussian(double* vars, double * data, int xlen, int ylen)
{
	double  tIstar = 0.;
	double* model;
	int j;
	GaussDataType * gdata;
	bfgs bfgs_o(target2DGaussian, 12); //optimisation object
	int osize = xlen * ylen;

	//reserve space for model
	if (model = (double*)malloc(sizeof(double) * osize)) {
		printf("error, out of memory\n");
		return -1;
	}

	//fill gdata struct
	gdata->data = data;
	gdata->model = model;
	gdata->xlen = xlen;
	gdata->ylen = ylen;

	//parameters 12-16 contain fit information and are fixed
	//1 Gauss fit uses first 6 parameters
	if ((int) vars[16] == 0) {
		for (j = 6; j < 12; j++) bfgs_o.fix(j);
	}
	//2 Gauss fit uses first 9 parameters
	else if ((int) vars[16] == 1) {
		for (j = 9; j < 12; j++) bfgs_o.fix(j);
	}
	//3 gauss fit uses first 12 parameters
	else if ((int) vars[16] == 2) {
		for (j = 12; j < 12; j++) bfgs_o.fix(j);
	}
	
	//set levenberg-marquadt conversion parametes
	//bfgs_o.seteps(0.001);
	bfgs_o.maxiter = 1000;

	//fix epsilon if indicated by function caller
	if (vars[15] == 1) { vars[4] = 1; bfgs_o.fix(4); }
	//fix bg if indicated by function caller
	if (vars[14] == 1) bfgs_o.fix(5);
	//if bg is free, make sure initial-guess is non-zero.
	//NV comment: do we really need this?
	if (vars[14] == 0 && vars[5] == 0) vars[5] = 0.1;

	vars[12] = bfgs_o.minimize(vars, gdata);

	//get magnitude of tIstar for optimised solution
	vars[17] = W2DG(gdata->data, gdata->model, osize);

	free(model);
	delete gdata;
	return 1;
}

/////////////////////////////////////////////////////////////////////////////////////////////

int Gauss2D_analysis_Ani(double* image,				// one frame
	int size_x, int size_y,		// image size
	OptsCluster* options,
	int Nimage,
	int& Nall, 			   //  total number of peaks found
	int& Ngood, 			//  total number of good peaks found
	int& Nconverged, 	//  total number of converged peaks
	int& Npeaks, 			//  number of remaining peaks 
	void** presults,		// see ResultsCluster definition
	int wi_nowi,			// weighted or no weights  ?
	int Npeaks_tmp,			// number of found peaks 
	int* peak_x_tmp,		// x coordinates of loaded peaks 
	int* peak_y_tmp,		// y coordinates of loaded peaks  
	int input_estimated_bg,	// with input or estimated bg value?
	MGParam* p)
{
	
	const int NPEAKS_MAX = options->maxNPeaks + 100;
	int offset;

	int i = 0, j, x, y;
	int yshift;
	
	double bg = 0., i0d;

	// first 4 bytes contain the array size -- must skip. See LabView help.
	ResultsCluster* results = (ResultsCluster*)((__int64*)(*presults) + 1);


	////////////////////////////////// search for close neighbours ////////////////////////////////

	int* badpeaks = new int[NPEAKS_MAX*NPEAKS_FACTOR];
	int ngoodpeaks = Npeaks_tmp;
	double distance;
	Nall += ngoodpeaks;

	for (i = 0; i < NPEAKS_MAX*NPEAKS_FACTOR; i++) badpeaks[i] = 0;

	for (i = 0; i < Npeaks_tmp; i++)
		for (j = i + 1; j < Npeaks_tmp; j++) {

			distance = sqrt((peak_x_tmp[i] - peak_x_tmp[j]) * (peak_x_tmp[i] - peak_x_tmp[j]) + (peak_y_tmp[i] - peak_y_tmp[j]) * (peak_y_tmp[i] - peak_y_tmp[j]));

			if (distance < (results[i].std + results[j].std + 1))
			{
				badpeaks[i] = 1;
				badpeaks[j] = 1;
				ngoodpeaks = ngoodpeaks - 2;
			}
		}

	if (ngoodpeaks > NPEAKS_MAX && ngoodpeaks <= 0) return 1;

	////////////////////////////////// centers of mass and fitting ////////////////////////////////////

	double cm_x, cm_y, sum_i, sigma_x, sigma_y;
	double* vars = new double[NVARS+4];

	// working arrays for LM fitting
//	double lm_options[] = { 1.e-10, 1.e-10, 1.e-10, 1.e-12, 200., 0.1, 1000.,  2. };  // default was { 1.e-10, 1.e-10, 1.e-10, 1.e-12, 100., 0.1, 200.,  2. }
	
	int info = 0;

	for (i = 0; i < Npeaks_tmp; i++) {

		offset = int(results[Npeaks].std);
		results[Npeaks].std = 0.;

		int osize = 2 * offset + 1;																		//object size
		int osize_sq = osize * osize;																	//object area
		double* icut = new double[osize_sq];

		double* f = new double[osize_sq];
//		double* weights = new double[osize_sq];

		if (badpeaks[i]) continue;
		else Ngood++;

		// copy
		j = 0;
		for (y = peak_y_tmp[i] - offset; y <= peak_y_tmp[i] + offset; y++) {
			yshift = size_x * y;
			for (x = peak_x_tmp[i] - offset; x <= peak_x_tmp[i] + offset; x++) {
				if ((0 <= (yshift + x)) && ((yshift + x) < size_x*size_y)) {
					icut[j] = image[yshift + x]; j++;
				}
				else {
					icut[j] = 0.; j++;
				}
			}
		}

//		save_array(icut, osize*osize);

		// bg estimate

			bg = 0.; int N_bg = 0;
			if (input_estimated_bg == 1) {
				x = peak_x_tmp[i] - offset - 1;
				for (y = peak_y_tmp[i] - offset; y <= peak_y_tmp[i] + offset; y++)
					if ((0 <= (size_x * y + x)) && ((size_x * y + x) < size_x*size_y)) { bg += image[size_x * y + x]; N_bg++; };		// left column

				y = peak_y_tmp[i] - offset - 1; yshift = size_x * y;
				for (x = peak_x_tmp[i] - offset - 1; x <= peak_x_tmp[i] + offset + 1; x++)
					if ((0 <= (yshift + x)) && ((yshift + x) < size_x*size_y)) { bg += image[yshift + x]; N_bg++; };					// top row

				x = peak_x_tmp[i] + offset + 1;
				for (y = peak_y_tmp[i] - offset; y <= peak_y_tmp[i] + offset; y++)
					if ((0 <= (size_x * y + x)) && ((size_x * y + x) < size_x*size_y)) { bg += image[size_x * y + x]; N_bg++; };		// right column

				y = peak_y_tmp[i] + offset + 1; yshift = size_x * y;
				for (x = peak_x_tmp[i] - offset - 1; x <= peak_x_tmp[i] + offset + 1; x++)
					if ((0 <= (yshift + x)) && ((yshift + x) < size_x*size_y)) { bg += image[yshift + x]; N_bg++; };					// bottom row
				bg /= N_bg;

			}
			else bg = options->background;

		// center of mass
		
		cm_x = 0.; cm_y = 0.; sum_i = 0.;
		for (y = 0; y < osize; y++) {
			yshift = osize * y;
			for (x = 0; x < osize; x++) {
				i0d = icut[yshift + x] - bg;
				cm_x += i0d * x;
				cm_y += i0d * y;
				sum_i += i0d;
			}
		};
		cm_x /= sum_i; cm_y /= sum_i;

				// estimat sigma of 2D Gauss
		double d_x, d_y;
		d_x = 0.; d_y = 0.;

		for (y = 0; y < osize; y++) {
			yshift = osize * y;
			for (x = 0; x < osize; x++) {
				i0d = icut[yshift + x];
				d_x += i0d * (x - cm_x)*(x - cm_x) / sum_i;
				d_y += i0d * (y - cm_y)*(y - cm_y) / sum_i;
			}
		}  

		results[Npeaks].sigma_x = sqrt(d_x);
		results[Npeaks].sigma_y = sqrt(d_y) / sqrt(d_x);


		// prepare to fit 2D Gauss
		// vars = [x0 y0 A sigma_x sigma_y bg info wi_nowi free_fixed_bg]: initial values
		//NV comment: this code needs to be changed to allow for fitting multiple Gaussians
		vars[0] = cm_x; 
		vars[1] = cm_y;
		vars[2] = icut[2 * offset*(offset + 1)] - bg;
		vars[3] = results[Npeaks].sigma_x; 
		vars[4] = results[Npeaks].sigma_y; 
		vars[5] = bg;
		vars[7] = (double) wi_nowi;
		vars[8] = (double)options->free_fixed;
		vars[9] = (double)options->elliptical_circular;

		// calculate weights
/*
		for (j = 0; j < osize_sq; j++)
		{
			if (wi_nowi == 0)
			{
				icut[j] <= 1. ? weights[j] = 1. : weights[j] = 1. / sqrt(icut[j]);
			}
			else weights[j] = 1.;
		}
*/
		/*
		comment NV for Suren: 
		I've changed the definition of fit2DGaussian.
		In the current definition there is no need to build the MGparam* struct.
		Therefore I've commented out the below piece of code and changed the call to 
		fit2DGaussian.
		Please verify that icut contains the data to be fitted.
		Also, fit2DGaussian now accepts non-square images.
		Please change the arguments (osize, osize) to (xlen, ylen)
		if you wish to use this function
		*/
		/*
		//copy data into LV cluster
		LVDoubleArray *subimage = *(p->subimage), *M = *(p->M);
		p->osize = osize;
		subimage->length = osize_sq;  
		M->length = osize_sq;

		for (j = 0; j < osize_sq; j++){
		subimage->data[j] = icut[j];
		M->data[j] = icut[j];
		}
		*/
		// fitting if requested
		if (options->fit2DGauss)
		{
			fit2DGaussian(vars, icut, osize, osize);
		}

		//NV comment: this code needs to be changed in order to accomodate multiple Gaussians
		// output
		for (int i = 0; i < 6; i++) vars[i] = fabs(vars[i]);
		sigma_x = fabs(vars[3]); sigma_y = fabs(vars[3] * vars[4]);

		// check sigma, intensity, and convergence
		if (options->fit2DGauss && options->must_converge && (vars[6] > 4)) continue;
		else Nconverged++;

//		if (fabs(sigma_x - sigma_y) / (sigma_x + sigma_y) > options->dsigma_max) continue;
//		if (options->relative_thr && (bg_surface != NULL) && (vars[2] < threshold)) continue;
//		if (fabs(vars[0]-offset)>offset || fabs(vars[1]-offset)>offset) continue;

		//NV comment: this code needs to be changed in order to accomodate multiple Gaussians
		results[Npeaks].imageID = Nimage;
		results[Npeaks].pixelID = results[Npeaks].max_y * size_x + results[Npeaks].max_x + 1;
		results[Npeaks].peak_x = peak_x_tmp[i] + vars[0] - offset;
		results[Npeaks].peak_y = peak_y_tmp[i] + vars[1] - offset;
		results[Npeaks].intensity = vars[2];
		results[Npeaks].sigma_x = sigma_x;
		results[Npeaks].sigma_y = sigma_y;
		results[Npeaks].background = vars[5];					// bg;
		results[Npeaks].max_x = peak_x_tmp[i];
		results[Npeaks].max_y = peak_y_tmp[i];
		
		sum_i = 0.;
		for (y = 0; y < osize; y++) {
			yshift = osize * y;
			for (x = 0; x < osize; x++) {
				sum_i += icut[yshift + x] - vars[5];	 
			}
		};

		results[Npeaks].Ncounts = sum_i;
		results[Npeaks].chi2 = osize_sq * vars[17];
		results[Npeaks].chi2 /= (osize_sq - NVARS + 1);
		results[Npeaks++].lm_message = (int)vars[6];

		delete[] f; 
//		delete[] weights;
		delete[] icut;
	}

	delete[] badpeaks; delete[] vars;

	return 0;
}


