// Fit 2DGauss
// Version 2019-04-16
//edited extensively by Nicolaas van der Voort
//AG Seidel, Düsseldorf

#include <iostream>
#include <fstream>
#include <math.h>
#include <set>
#include "fit2DGaussian.h"
#include "i_lbfgs.h"
#include "twoIstar.h"

using namespace std;

const int NVARS = 6;
const int NPEAKS_FACTOR = 10;

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


///////////////////////////////// target function (to minimize) ///////////////////////////////

double target2DGaussian(double* vars, void* pM)
{
	double w;
	MGParam* p = (MGParam*)pM;
	LVDoubleArray *subimage = *(p->subimage), *M = *(p->M);
	std::set<std::int8_t> positionSet = { 0, 1, 6, 7, 9, 10 };

	//check that the initial positions of all gaussians remain within the frame
	std::set<std::int8_t>::iterator it = positionSet.begin();
	while (it != positionSet.end()) {
		if ((vars[*it] < 0.) || (vars[*it] > p->osize - 1)) vars[*it] = (p->osize) / 2.;
		it++;
	}
 
	//get model
	if ((int) vars[16] == 0 ) {
		model2DGaussian(vars, M->data, p->osize);
	}
	else if ((int) vars[16] == 1) {
		modelTwo2DGaussian(vars, M->data, p->osize);
	}
	else if ((int) vars[16] == 2) {
		modelThree2DGaussian(vars, M->data, p->osize);
	}

	w = W2DG(subimage->data, M->data, (int)vars[7], subimage->length);
	return w;
}

//////////////////////////////////////////// model2DGaussian ////////////////////////////////////////////
// vars = [x0 y0 A sigma ellipticity bg]
int model2DGaussian(double* vars, double* M, int osize)
{
	int x, y, j = 0;
	for (int i = 0; i < 6; i++) vars[i] = fabs(vars[i]);
	double x0 = vars[0], y0 = vars[1], A = vars[2], B = vars[5];
	double tx = 0.5 / (vars[3] * vars[3]), ty = 0.5 / (vars[3] * vars[3] * vars[4] * vars[4]);
	double* ex = new double[osize]; double ey;

	for (x = 0; x < osize; x++) ex[x] = exp(-(x - x0)*(x - x0)*tx);

	for (y = 0; y < osize; y++) {
		ey = exp(-(y - y0)*(y - y0)*ty);
		for (x = 0; x < osize; x++)
		{
			M[j] = (A*ex[x] * ey + B); j++;					// Model of 2D Gaussian
		}
	}

	delete[] ex;
	return 0;
}

//////////////////////////////////////////// modelTwo2DGaussian ////////////////////////////////////////////
//function uses model2DGaussian function for constructor, see fit2DGaussian for parameter declaration
int modelTwo2DGaussian(double* vars, double* M, int osize) {
	double * vars_dummy = new double[6];
	double * M_dummy = new double[osize * osize];
	int i, osize_sq = osize * osize;
	//create the first Gaussian with bg, store in M
	model2DGaussian(vars, M, osize);

	//create the second Gaussian without bg, store in M_dummy
	vars_dummy[0] = vars[6];
	vars_dummy[1] = vars[7];
	vars_dummy[2] = vars[8];
	vars_dummy[3] = vars[3];
	vars_dummy[4] = vars[4];
	vars_dummy[5] = 0;
	model2DGaussian(vars_dummy, M_dummy, osize);

	//add
	for (i = 0; i < osize_sq; i++) {
		M[i] += M_dummy[i];
	}

	delete[] vars_dummy, M_dummy;
	return 0;
}

//////////////////////////////////////////// modelThree2DGaussian ////////////////////////////////////////////

//function uses model2DGaussian and model2DGaussian function for constructor, see fit2DGaussian for parameter declaration
int modelThree2DGaussian(double* vars, double* M, int osize) {
	double * vars_dummy = new double[6];
	double * M_dummy = new double[osize * osize];
	int i, osize_sq = osize * osize;
	//create the first two Gaussian with bg, store in M
	modelTwo2DGaussian(vars, M, osize);

	//create the third Gaussian without bg, store in M_dummy
	vars_dummy[0] = vars[9];
	vars_dummy[1] = vars[10];
	vars_dummy[2] = vars[11];
	vars_dummy[3] = vars[3];
	vars_dummy[4] = vars[4];
	vars_dummy[5] = 0;

	model2DGaussian(vars_dummy, M_dummy, osize);

	//add
	for (i = 0; i < osize_sq; i++) {
		M[i] += M_dummy[i];
	}

	delete[] vars_dummy, M_dummy;
	return 0;
}

//////////////////////////////////////////// fit2DGaussian ////////////////////////////////////////////

//fit2DGaussian initializes optimisation routine
//vars contains the parameters that are optimized
//  vars = [0: x0, 1: y0, 2: A0, 3: sigma, 4: ellipticity, 5: bg, 6: x1, 7: y1, 8: A1, 9: x2, 10: y2, 11: A2, 12: info,\
//          13: wi_nowi, 14: fit_bg, 15: ellipt_circ, 16: model]
//  info contains information from the fitting algorithm
//  wi_nowi contains weights or no_wights, outdated?
//  fit_bg asks if background is fitted. 0 -> bg is fitted
//  ellipt_circ  determines if elliptical fits are allowed. Only relevant for model2DGaussian
//  model determines the model to be used:
//    0: model2DGaussian
//    1: modelTwo2DGaussian
//    2: modelThree2DGaussian
//MGParam Class contains the data function, the model function and the number of rows in the image.
/*
Eigen::VectroXd fitNew(const Eigen::MatrixXi &img)
{
	Eigen::VectroXd fit;
	
	return fit;
}
*/
double fit2DGaussian(double* vars, MGParam* p)
{
	double  tIstar = 0.;
	LVDoubleArray *subimage = *(p->subimage), *M = *(p->M);
	int j;

	//fix parameters according to model
	bfgs bfgs_o(target2DGaussian, 17);
	//bfgs_o.seteps(0.001);
	if ((int) vars[16] == 0) {
		for (j = 6; j < 17; j++) bfgs_o.fix(j);
	}
	else if ((int) vars[16] == 1) {
		for (j = 9; j < 17; j++) bfgs_o.fix(j);
	}
	else if ((int) vars[16] == 2) {
		for (j = 12; j < 17; j++) bfgs_o.fix(j);
	} 

	bfgs_o.maxiter = 1000;
	if (vars[15] == 1) { vars[4] = 1; bfgs_o.fix(4); }
	if (vars[14] == 1) bfgs_o.fix(5);
	//NV comment: do we really need this?
	if (vars[14] == 0 && vars[5] == 0) vars[5] = 0.1;
	vars[12] = bfgs_o.minimize(vars, p);

	//get magnitude of tIstar for optimised solution
	tIstar = W2DG(subimage->data, M->data, (int)vars[7], subimage->length);
	return tIstar;
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
	
	double bg = 0., i0d, tIstar;

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
		//copy data into LV cluster
		LVDoubleArray *subimage = *(p->subimage), *M = *(p->M);
		p->osize = osize;
		subimage->length = osize_sq;  
		M->length = osize_sq;

		for (j = 0; j < osize_sq; j++){
		subimage->data[j] = icut[j];
		M->data[j] = icut[j];
		}

		// fitting if requested
		if (options->fit2DGauss)
		{
			tIstar = fit2DGaussian(vars, p);
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
		results[Npeaks].chi2 = osize_sq * tIstar;
		results[Npeaks].chi2 /= (osize_sq - NVARS + 1);
		results[Npeaks++].lm_message = (int)vars[6];

		delete[] f; 
//		delete[] weights;
		delete[] icut;
	}

	delete[] badpeaks; delete[] vars;

	return 0;
}


