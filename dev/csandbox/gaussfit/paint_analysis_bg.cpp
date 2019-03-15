// paint analysis v. 2008-11-14
#ifdef BACKGROUND_FIXED

#include <math.h>
#include "paint.h"

const int NVARS = 5;
const int NPEAKS_FACTOR = 10;

// background and threshold surfaces
static double* bg_surface = NULL;
static int* threshold_surface = NULL;
static int last_threshold = -108;

// weighted residuals, weights*(model-data)

void wr_2DGauss(double* f, double* vars, int osize, double* icut, double* weights, double ABscale, double bg)
{
  // vars = [x0 y0 A sigma_x sigma_y]

  int x,y,j=0;
  double x0 = vars[0], y0 = vars[1], A = vars[2]*ABscale, B = bg*ABscale;
  double tx = 0.5/(vars[3]*vars[3]), ty = 0.5/(vars[4]*vars[4]);
  double* ex = new double[osize]; double ey;
  for (x = 0; x < osize; x++) ex[x] = exp(-(x-x0)*(x-x0)*tx);

  for (y = 0; y < osize; y++) {
    ey = exp(-(y-y0)*(y-y0)*ty);
    for (x = 0; x < osize; x++)
      f[j++] = (A*ex[x]*ey+B-icut[j])*weights[j];
  }
  delete [] ex;
}

// gradient of wr_2DGauss

void wr_2DGauss_grad(double* fjac, double* vars, int osize, double* weights, double ABscale, double bg)
{
  // vars = [x0 y0 A sigma_x sigma_y]

  int x,y,j=0;
  double x0 = vars[0], y0 = vars[1];
  double A = vars[2]*ABscale, ew, tx = 0.5/(vars[3]*vars[3]), ty = 0.5/(vars[4]*vars[4]),
    dsx = tx*2./vars[3], dsy = ty*2./vars[4];
  int osize_sq = osize*osize;
  double* ex = new double[osize]; double ey;
  for (x = 0; x < osize; x++) ex[x] = exp(-(x-x0)*(x-x0)*tx);

  for (y = 0; y < osize; y++) {
    ey = exp(-(y-y0)*(y-y0)*ty);
    for (x = 0; x < osize; x++) {
      j = y*osize+x;     
      ew = ex[x]*ey*weights[j];

      fjac[j] = A*ew*tx*2.*(x-x0);			// d/dx0
      fjac[j+osize_sq] = A*ew*ty*2.*(y-y0);		// d/dy0
      fjac[j+2*osize_sq] = ABscale*ew;			// d/dA
      fjac[j+3*osize_sq] = A*ew*dsx*(x-x0)*(x-x0);	// d/dsigma_x
      fjac[j+4*osize_sq] = A*ew*dsy*(y-y0)*(y-y0);	// d/dsigma_y
    }
  }
  delete [] ex;
}


/////////////////////////////////////////////////////////////////////////////////////////////

int paint_analysis(int* image,				// one frame
			int size_x, int size_y,		// image size
			OptsCluster* options,
			int Nimage,
			int& Nall, 			   //  total number of peaks found
			int& Ngood, 			//  total number of good peaks found
			int& Nconverged, 	//  total number of converged peaks
			int& Npeaks, 			//  number of remaining peaks 
			void** presults,		// see ResultsCluster definition
			int wi_nowi)      // weighted or no weights  ?
{

 const int NPEAKS_MAX = options->maxNPeaks;
 const int offset = options->offset;

 int i=0, j, x, y;
 int yshift, xy;
 int i0;
 double bg = 0., i0d, ABscale;

 // first 4 bytes contain the array size -- must skip. See LabView help.
 ResultsCluster* results = (ResultsCluster*)((__int64*)(*presults)+1);

 //////////////////////////////////////// find peaks ///////////////////////////////////////////

 int Npeaks_tmp = 0;
 int* peak_x_tmp = new int[NPEAKS_MAX*NPEAKS_FACTOR];
 int* peak_y_tmp = new int[NPEAKS_MAX*NPEAKS_FACTOR];


 // init threshold_surface
 const int threshold = options->threshold;
 const int max_threshold = options->max_threshold;
 if (options->relative_thr && (threshold!=last_threshold) && (bg_surface!=NULL)) {
   for(y = 0; y < size_y; y++)
     for(x = 0; x < size_x; x++)
       threshold_surface[size_x*y+x] = threshold + (int)bg_surface[size_x*y+x];
   last_threshold = threshold; 
 }
 
  if (options->relative_thr && (bg_surface!=NULL)) // search with respect to bg_surface

   for(y = offset; y < size_y-offset; y++) {
     yshift = size_x*y;
     for(x = offset; x < size_x-offset; x++) {
       xy = yshift+x;
       i0 = image[xy];
       if ((i0 >= threshold_surface[xy]) && 
		   ((image[xy - 1]< i0) & (image[xy + 1] <= i0)) &&
		   ((image[xy - size_x - 1] < i0) & (image[xy - size_x] < i0) & (image[xy - size_x + 1] < i0)) &&
		   ((image[xy + size_x - 1] < i0) & (image[xy + size_x] <= i0) & (image[xy + size_x + 1] <= i0)))
       {
         peak_x_tmp[Npeaks_tmp] = x;
         peak_y_tmp[Npeaks_tmp++] = y;
         if (Npeaks_tmp >= NPEAKS_MAX*NPEAKS_FACTOR) return 1;
       }
     }
   
   }


 else 		// normal search, value > threshold

   for(y = offset; y < size_y-offset; y++) {
     yshift = size_x*y;
     for(x = offset; x < size_x-offset; x++) {
       xy = yshift+x;
       i0 = image[xy];
       if ((i0 >= threshold) && (i0 <= max_threshold) &&
          ((image[xy-1]< i0) & (image[xy+1] <= i0)) &&
          ((image[xy-size_x-1] < i0) & (image[xy-size_x] < i0) & (image[xy-size_x+1] < i0)) &&
          ((image[xy+size_x-1] < i0) & (image[xy+size_x] <= i0) & (image[xy+size_x+1] <= i0)))
       {
         peak_x_tmp[Npeaks_tmp] = x;
         peak_y_tmp[Npeaks_tmp++] = y;
         if (Npeaks_tmp >= NPEAKS_MAX*NPEAKS_FACTOR) return 1;
       }
     }
     
   }
   
 ////////////////////////////////// search for close neighbours ////////////////////////////////

 int* badpeaks = new int[NPEAKS_MAX*NPEAKS_FACTOR];
 int ngoodpeaks = Npeaks_tmp;
 Nall += ngoodpeaks;
 
 ///////
 
 for (i = 0; i < NPEAKS_MAX*NPEAKS_FACTOR; i++) badpeaks[i] = 0;

 for (i = 0; i < Npeaks_tmp; i++) 
   for (j = i+1; j < Npeaks_tmp; j++)
     if ((abs((int)peak_x_tmp[i]-(int)peak_x_tmp[j]) < 2*offset) &
         (abs((int)peak_y_tmp[i]-(int)peak_y_tmp[j]) < 2*offset)) {
       badpeaks[i] = 1;
       badpeaks[j] = 1;
       ngoodpeaks = ngoodpeaks - 2;
     }
 if (ngoodpeaks > NPEAKS_MAX && ngoodpeaks <=0) return 1;
 
 ////////////////////////////////// centers of mass and fitting ////////////////////////////////////
 
 int osize = 2*offset + 1;
 int osize_sq = osize*osize;
 double* icut = new double[osize_sq];
 
 double cm_x, cm_y, sum_i, sigma_x, sigma_y;

 // working arrays for LM fitting
 
 double* vars = new double[NVARS];
 double* f = new double[osize_sq];
 double* weights = new double[osize_sq];
 int lmresult = 0;

 double lm_options[] = {1.e-10, 1.e-10, 1.e-10, 1.e-12, 100., 0.1, 200.,  2.};
 double *diag, *fjac, *qtf, *wa1, *wa2, *wa3, *wa4;
 int *ipvt;
 diag = new double[NVARS];
 for (j = 0; j < NVARS; j++) diag[j] = 1.;
 fjac = new double[NVARS*osize_sq];
 qtf = new double[NVARS];
 wa1 = new double[NVARS];
 wa2 = new double[NVARS];
 wa3 = new double[NVARS];
 wa4 = new double[osize_sq];
 ipvt  = new int[NVARS];

 int* s_int = new int[8];
 double* s_dbl = new double[16];

 int iflag=0, nfev=0, info=10;
 
 for (i = 0; i < Npeaks_tmp; i++) {

   if (badpeaks[i]) continue;
   else Ngood++;
   
   // copy
   j = 0;
   for (y = peak_y_tmp[i]-offset; y <= peak_y_tmp[i]+offset; y++) {
     yshift = size_x*y;
     for (x = peak_x_tmp[i]-offset; x <= peak_x_tmp[i]+offset; x++)
       icut[j++] = image[yshift+x];
   }

   // bg estimate
   if (options->relative_thr && (bg_surface!=NULL))
     bg = bg_surface[size_x*peak_y_tmp[i] + peak_x_tmp[i]];
   else {
     bg = 0.;
     for (x = 0; x < osize; x++) bg += icut[x];
     for (x = osize_sq - osize; x < osize_sq; x++) bg += icut[x];
     bg /= (2.*(double)osize);
   }

   // center of mass
   cm_x = 0.; cm_y = 0.; sum_i = 0.;
   for (y = 0; y < osize; y++) {
     yshift = osize*y;
     for (x = 0; x < osize; x++) {
       i0d = icut[yshift+x] - bg;
       cm_x += i0d * x;
       cm_y += i0d * y;
       sum_i += i0d;
     }
   }
   cm_x /= sum_i; cm_y /= sum_i;

   
   // prepare to fit 2D Gauss

   // vars = [x0 y0 A B sigma_x sigma_y]: initial values
   double bg_fit;
   vars[0] = cm_x; vars[1] = cm_y; 
   ABscale = icut[2*offset*(offset+1)] - bg; vars[2] = 1.; bg_fit = bg/ABscale;
   vars[3] = options->sigma0; vars[4] = options->sigma0;

   // calculate weights and eval f
   
   for (j = 0; j < osize_sq; j++)
   {
	   if (wi_nowi==0)
	   {
      icut[j]<=1. ? weights[j] = 1. : weights[j] = 1./sqrt(icut[j]);
	   }
	   else weights[j] = 1.; 
   }
   wr_2DGauss(f, vars, osize, icut, weights, ABscale, bg_fit);

   // fitting if requested

   if (options->fit2DGauss) {

     // init lmdif
     for (j = 0; j < 8; j++) s_int[j] = 0;
     for (j = 0; j < 16; j++) s_dbl[j] = 0.;
     nfev = 0; info = 0;
 
     do {
       lmresult = lmdif(vars, f, lm_options,
            osize_sq, NVARS, osize_sq,
            diag, fjac, ipvt, qtf,
            wa1, wa2, wa3, wa4,
            s_int, s_dbl, &info, &nfev);
       if (lmresult==1) wr_2DGauss(f, vars, osize, icut, weights, ABscale, bg_fit);
       if (lmresult==2) {
         wr_2DGauss(f, vars, osize, icut, weights, ABscale, bg_fit);
         wr_2DGauss_grad(fjac, vars, osize, weights, ABscale, bg_fit);
       }
     }
     while (lmresult);
   }
     
   // output

   // check sigma, intensity, and convergence
  
   sigma_x = fabs(vars[3]); sigma_y = fabs(vars[4]); 
   
   
   if (options->fit2DGauss && options->must_converge && (info >= 4)) continue;
   else Nconverged++;
   
   if (fabs(sigma_x-sigma_y)/(sigma_x+sigma_y) > options->dsigma_max) continue;  

   if (options->relative_thr && (bg_surface!=NULL) && (vars[2]*ABscale < threshold)) continue;
   // if (abs(vars[0]-offset)>offset || abs(vars[1]-offset)>offset) continue;
   results[Npeaks].imageID = Nimage;
   results[Npeaks].pixelID = peak_y_tmp[i]*size_x + peak_x_tmp[i];
   results[Npeaks].peak_x = peak_x_tmp[i] + vars[0] - offset;
   results[Npeaks].peak_y = peak_y_tmp[i] + vars[1] - offset;
   results[Npeaks].intensity = vars[2]*ABscale;
   results[Npeaks].sigma_x = sigma_x;
   results[Npeaks].sigma_y = sigma_y;
   results[Npeaks].background = bg_fit*ABscale;
   results[Npeaks].max_x = peak_x_tmp[i];
   results[Npeaks].max_y = peak_y_tmp[i];   
   results[Npeaks].Ncounts = sum_i; 
   results[Npeaks].chi2 = 0.;
   for (j = 0; j < osize_sq; j++) results[Npeaks].chi2 += f[j]*f[j];
   results[Npeaks].chi2 /= (osize_sq - NVARS + 1);
   results[Npeaks++].lm_message = info;
 }


 delete [] peak_x_tmp; delete [] peak_y_tmp; delete [] badpeaks;
 delete [] vars; delete [] f; delete [] weights;
 delete [] icut; delete [] diag; delete [] fjac; delete [] qtf; 
 delete [] wa1; delete [] wa2; delete [] wa3; delete [] wa4; 
 delete [] ipvt; delete [] s_int; delete [] s_dbl;

 return 0;

}

/////////////////////////////////////////////////////////////////////////////////////////////

int Gauss2D_analysis(int* image,				// one frame
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
	int* peak_y_tmp)     // y coordinates of loaded peaks  
{

	const int NPEAKS_MAX = options->maxNPeaks;
	const int offset = options->offset;

	int i = 0, j, x, y;
	int yshift, xy;
	int i0;
	double bg = 0., i0d, ABscale;

	// first 4 bytes contain the array size -- must skip. See LabView help.
	ResultsCluster* results = (ResultsCluster*)((__int64*)(*presults) + 1);

	//////////////////////////////////////// find peaks ///////////////////////////////////////////
/*
	int Npeaks_tmp = 0;
	int* peak_x_tmp = new int[NPEAKS_MAX*NPEAKS_FACTOR];
	int* peak_y_tmp = new int[NPEAKS_MAX*NPEAKS_FACTOR];

*/
	// init threshold_surface
	const int threshold = options->threshold;
	const int max_threshold = options->max_threshold;
	if (options->relative_thr && (threshold != last_threshold) && (bg_surface != NULL)) {
		for (y = 0; y < size_y; y++)
			for (x = 0; x < size_x; x++)
				threshold_surface[size_x*y + x] = threshold + (int)bg_surface[size_x*y + x];
		last_threshold = threshold;
	}
/*
	if (options->relative_thr && (bg_surface != NULL)) // search with respect to bg_surface

		for (y = offset; y < size_y - offset; y++) {
			yshift = size_x * y;
			for (x = offset; x < size_x - offset; x++) {
				xy = yshift + x;
				i0 = image[xy];
				if ((i0 >= threshold_surface[xy]) &&
					((image[xy - 1]< i0) & (image[xy + 1] <= i0)) &&
					((image[xy - size_x - 1] < i0) & (image[xy - size_x] < i0) & (image[xy - size_x + 1] < i0)) &&
					((image[xy + size_x - 1] < i0) & (image[xy + size_x] <= i0) & (image[xy + size_x + 1] <= i0)))
				{
					peak_x_tmp[Npeaks_tmp] = x;
					peak_y_tmp[Npeaks_tmp++] = y;
					if (Npeaks_tmp >= NPEAKS_MAX * NPEAKS_FACTOR) return 1;
				}
			}

		}


	else 		// normal search, value > threshold

		for (y = offset; y < size_y - offset; y++) {
			yshift = size_x * y;
			for (x = offset; x < size_x - offset; x++) {
				xy = yshift + x;
				i0 = image[xy];
				if ((i0 >= threshold) && (i0 <= max_threshold) &&
					((image[xy - 1]< i0) & (image[xy + 1] <= i0)) &&
					((image[xy - size_x - 1] < i0) & (image[xy - size_x] < i0) & (image[xy - size_x + 1] < i0)) &&
					((image[xy + size_x - 1] < i0) & (image[xy + size_x] <= i0) & (image[xy + size_x + 1] <= i0)))
				{
					peak_x_tmp[Npeaks_tmp] = x;
					peak_y_tmp[Npeaks_tmp++] = y;
					if (Npeaks_tmp >= NPEAKS_MAX * NPEAKS_FACTOR) return 1;
				}
			}

		}
		*/

	////////////////////////////////// search for close neighbours ////////////////////////////////

	int* badpeaks = new int[NPEAKS_MAX*NPEAKS_FACTOR];
	int ngoodpeaks = Npeaks_tmp;
	Nall += ngoodpeaks;

	///////

	for (i = 0; i < NPEAKS_MAX*NPEAKS_FACTOR; i++) badpeaks[i] = 0;

	for (i = 0; i < Npeaks_tmp; i++)
		for (j = i + 1; j < Npeaks_tmp; j++)
			if ((abs((int)peak_x_tmp[i] - (int)peak_x_tmp[j]) < 2 * offset) &
				(abs((int)peak_y_tmp[i] - (int)peak_y_tmp[j]) < 2 * offset)) {
				badpeaks[i] = 1;
				badpeaks[j] = 1;
				ngoodpeaks = ngoodpeaks - 2;
			}
	if (ngoodpeaks > NPEAKS_MAX && ngoodpeaks <= 0) return 1;

	////////////////////////////////// centers of mass and fitting ////////////////////////////////////

	int osize = 2 * offset + 1;
	int osize_sq = osize * osize;
	double* icut = new double[osize_sq];

	double cm_x, cm_y, sum_i, sigma_x, sigma_y;

	// working arrays for LM fitting

	double* vars = new double[NVARS];
	double* f = new double[osize_sq];
	double* weights = new double[osize_sq];
	int lmresult = 0;

	double lm_options[] = { 1.e-10, 1.e-10, 1.e-10, 1.e-12, 100., 0.1, 200.,  2. };
	double *diag, *fjac, *qtf, *wa1, *wa2, *wa3, *wa4;
	int *ipvt;
	diag = new double[NVARS];
	for (j = 0; j < NVARS; j++) diag[j] = 1.;
	fjac = new double[NVARS*osize_sq];
	qtf = new double[NVARS];
	wa1 = new double[NVARS];
	wa2 = new double[NVARS];
	wa3 = new double[NVARS];
	wa4 = new double[osize_sq];
	ipvt = new int[NVARS];

	int* s_int = new int[8];
	double* s_dbl = new double[16];

	int iflag = 0, nfev = 0, info = 10;

	for (i = 0; i < Npeaks_tmp; i++) {

		if (badpeaks[i]) continue;
		else Ngood++;

		// copy
		j = 0;
		for (y = peak_y_tmp[i] - offset; y <= peak_y_tmp[i] + offset; y++) {
			yshift = size_x * y;
			for (x = peak_x_tmp[i] - offset; x <= peak_x_tmp[i] + offset; x++)
				icut[j++] = image[yshift + x];
		}

		// bg estimate
		if (options->relative_thr && (bg_surface != NULL))
			bg = bg_surface[size_x*peak_y_tmp[i] + peak_x_tmp[i]];
		else {
			bg = 0.;
			for (x = 0; x < osize; x++) bg += icut[x];
			for (x = osize_sq - osize; x < osize_sq; x++) bg += icut[x];
			bg /= (2.*(double)osize);
		}

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
		}
		cm_x /= sum_i; cm_y /= sum_i;


		// prepare to fit 2D Gauss

		// vars = [x0 y0 A B sigma_x sigma_y]: initial values
		double bg_fit;
		vars[0] = cm_x; vars[1] = cm_y;
		ABscale = icut[2 * offset*(offset + 1)] - bg; vars[2] = 1.; bg_fit = bg / ABscale;
		vars[3] = options->sigma0; vars[4] = options->sigma0;

		// calculate weights and eval f

		for (j = 0; j < osize_sq; j++)
		{
			if (wi_nowi == 0)
			{
				icut[j] <= 1. ? weights[j] = 1. : weights[j] = 1. / sqrt(icut[j]);
			}
			else weights[j] = 1.;
		}
		wr_2DGauss(f, vars, osize, icut, weights, ABscale, bg_fit);

		// fitting if requested

		if (options->fit2DGauss) 
		{

			// init lmdif
			for (j = 0; j < 8; j++) s_int[j] = 0;
			for (j = 0; j < 16; j++) s_dbl[j] = 0.;
			nfev = 0; info = 0;

			do {
				lmresult = lmdif(vars, f, lm_options,
					osize_sq, NVARS, osize_sq,
					diag, fjac, ipvt, qtf,
					wa1, wa2, wa3, wa4,
					s_int, s_dbl, &info, &nfev);
				if (lmresult == 1) wr_2DGauss(f, vars, osize, icut, weights, ABscale, bg_fit);
				if (lmresult == 2) {
					wr_2DGauss(f, vars, osize, icut, weights, ABscale, bg_fit);
					wr_2DGauss_grad(fjac, vars, osize, weights, ABscale, bg_fit);
				}
			} while (lmresult);
		}

		// output

		// check sigma, intensity, and convergence

		sigma_x = fabs(vars[3]); sigma_y = fabs(vars[4]);


		if (options->fit2DGauss && options->must_converge && (info >= 4)) continue;
		else Nconverged++;

		if (fabs(sigma_x - sigma_y) / (sigma_x + sigma_y) > options->dsigma_max) continue;

		if (options->relative_thr && (bg_surface != NULL) && (vars[2] * ABscale < threshold)) continue;
		// if (abs(vars[0]-offset)>offset || abs(vars[1]-offset)>offset) continue;
		results[Npeaks].imageID = Nimage;
		results[Npeaks].pixelID = peak_y_tmp[i] * size_x + peak_x_tmp[i];
		results[Npeaks].peak_x = peak_x_tmp[i] + vars[0] - offset;
		results[Npeaks].peak_y = peak_y_tmp[i] + vars[1] - offset;
		results[Npeaks].intensity = vars[2] * ABscale;
		results[Npeaks].sigma_x = sigma_x; 
		results[Npeaks].sigma_y = sigma_y;
		results[Npeaks].background = bg_fit * ABscale;
		results[Npeaks].max_x = peak_x_tmp[i];
		results[Npeaks].max_y = peak_y_tmp[i];
		results[Npeaks].Ncounts = sum_i;
		results[Npeaks].chi2 = 0.;
		for (j = 0; j < osize_sq; j++) results[Npeaks].chi2 += f[j] * f[j];
		results[Npeaks].chi2 /= (osize_sq - NVARS + 1);
		results[Npeaks++].lm_message = info;
	}


//	delete[] peak_x_tmp; delete[] peak_y_tmp; 
	delete[] badpeaks;
	delete[] vars; delete[] f; delete[] weights;
	delete[] icut; delete[] diag; delete[] fjac; delete[] qtf;
	delete[] wa1; delete[] wa2; delete[] wa3; delete[] wa4;
	delete[] ipvt; delete[] s_int; delete[] s_dbl;

	return 0;

}

/////////////////////////////////////////////////////////////////////////////////////////////

int Gauss2D_analysis2(int* image,				// one frame
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
	int* peak_y_tmp)     // y coordinates of loaded peaks  
{

	const int NPEAKS_MAX = options->maxNPeaks+100;
	int offset = options->offset;

	int i = 0, j, x, y;
	int yshift, xy;
	int i0;
	double bg = 0., i0d, ABscale;

	// first 4 bytes contain the array size -- must skip. See LabView help.
	ResultsCluster* results = (ResultsCluster*)((__int64*)(*presults) + 1);

	//////////////////////////////////////// init threshold_surface ///////////////////////////////////////////

	const int threshold = options->threshold;
	const int max_threshold = options->max_threshold;
	if (options->relative_thr && (threshold != last_threshold) && (bg_surface != NULL)) {
		for (y = 0; y < size_y; y++)
			for (x = 0; x < size_x; x++)
				threshold_surface[size_x*y + x] = threshold + (int)bg_surface[size_x*y + x];
		last_threshold = threshold;
	}
	
	////////////////////////////////// search for close neighbours ////////////////////////////////

	int* badpeaks = new int[NPEAKS_MAX*NPEAKS_FACTOR];
	int ngoodpeaks = Npeaks_tmp;
	Nall += ngoodpeaks;

	for (i = 0; i < NPEAKS_MAX*NPEAKS_FACTOR; i++) badpeaks[i] = 0;

	for (i = 0; i < Npeaks_tmp; i++)
		for (j = i + 1; j < Npeaks_tmp; j++)
			if ((abs((int)peak_x_tmp[i] - (int)peak_x_tmp[j]) < 2 * offset) &
				(abs((int)peak_y_tmp[i] - (int)peak_y_tmp[j]) < 2 * offset)) {
				badpeaks[i] = 1;
				badpeaks[j] = 1;
				ngoodpeaks = ngoodpeaks - 2;
			}
	if (ngoodpeaks > NPEAKS_MAX && ngoodpeaks <= 0) return 1;

	////////////////////////////////// centers of mass and fitting ////////////////////////////////////


	double cm_x, cm_y, sum_i, sigma_x, sigma_y;

	// working arrays for LM fitting

	double* vars = new double[NVARS];

	int lmresult = 0;

	double lm_options[] = { 1.e-10, 1.e-10, 1.e-10, 1.e-12, 100., 0.1, 200.,  2. };
	double *diag, *fjac, *qtf, *wa1, *wa2, *wa3, *wa4;
	int *ipvt;
	diag = new double[NVARS];
	for (j = 0; j < NVARS; j++) diag[j] = 1.;

	qtf = new double[NVARS];
	wa1 = new double[NVARS];
	wa2 = new double[NVARS];
	wa3 = new double[NVARS];
	
	ipvt = new int[NVARS];

	int* s_int = new int[8];
	double* s_dbl = new double[16];

	int iflag = 0, nfev = 0, info = 10;

	for (i = 0; i < Npeaks_tmp; i++) {
		if (results[Npeaks].sigma_x >= results[Npeaks].sigma_y) offset = results[Npeaks].sigma_x;
		else offset = results[Npeaks].sigma_y;
		
		int osize = 2 * offset + 1;
		int osize_sq = osize * osize;
		double* icut = new double[osize_sq];
		
		double* f = new double[osize_sq];
		double* weights = new double[osize_sq];
		fjac = new double[NVARS*osize_sq];
		wa4 = new double[osize_sq];


		if (badpeaks[i]) continue;
		else Ngood++;

		// copy
		j = 0;
		for (y = peak_y_tmp[i] - offset; y <= peak_y_tmp[i] + offset; y++) {
			yshift = size_x * y;
			for (x = peak_x_tmp[i] - offset; x <= peak_x_tmp[i] + offset; x++)
				if ((0 <= (yshift + x)) && ((yshift + x) < size_x*size_y)) icut[j++] = image[yshift + x];
				else icut[j++] = 0;
		}

		// bg estimate
		if (options->relative_thr && (bg_surface != NULL))
			bg = bg_surface[size_x*peak_y_tmp[i] + peak_x_tmp[i]];
		else {
			bg = 0.;
			for (x = 0; x < osize; x++) bg += icut[x];
			for (x = osize_sq - osize; x < osize_sq; x++) bg += icut[x];
			bg /= (2.*(double)osize);
		}

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
		}
		cm_x /= sum_i; cm_y /= sum_i;


		// prepare to fit 2D Gauss

		// vars = [x0 y0 A B sigma_x sigma_y]: initial values
		double bg_fit;
		vars[0] = cm_x; vars[1] = cm_y;
		ABscale = icut[2 * offset*(offset + 1)] - bg; vars[2] = 1.; bg_fit = bg / ABscale;
		vars[3] = results[Npeaks].sigma_x; vars[4] = results[Npeaks].sigma_y;

		// calculate weights and eval f

		for (j = 0; j < osize_sq; j++)
		{
			if (wi_nowi == 0)
			{
				icut[j] <= 1. ? weights[j] = 1. : weights[j] = 1. / sqrt(icut[j]);
			}
			else weights[j] = 1.;
		}
		wr_2DGauss(f, vars, osize, icut, weights, ABscale, bg_fit);

		// fitting if requested

		if (options->fit2DGauss)
		{

			// init lmdif
			for (j = 0; j < 8; j++) s_int[j] = 0;
			for (j = 0; j < 16; j++) s_dbl[j] = 0.;
			nfev = 0; info = 0;

			do {
				lmresult = lmdif(vars, f, lm_options,
					osize_sq, NVARS, osize_sq,
					diag, fjac, ipvt, qtf,
					wa1, wa2, wa3, wa4,
					s_int, s_dbl, &info, &nfev);
				if (lmresult == 1) wr_2DGauss(f, vars, osize, icut, weights, ABscale, bg_fit);
				if (lmresult == 2) {
					wr_2DGauss(f, vars, osize, icut, weights, ABscale, bg_fit);
					wr_2DGauss_grad(fjac, vars, osize, weights, ABscale, bg_fit);
				}
			} while (lmresult);
		}

		// output

		// check sigma, intensity, and convergence

		sigma_x = fabs(vars[3]); sigma_y = fabs(vars[4]);


		if (options->fit2DGauss && options->must_converge && (info >= 4)) continue;
		else Nconverged++;

		if (fabs(sigma_x - sigma_y) / (sigma_x + sigma_y) > options->dsigma_max) continue;

		if (options->relative_thr && (bg_surface != NULL) && (vars[2] * ABscale < threshold)) continue;

		// if (abs(vars[0]-offset)>offset || abs(vars[1]-offset)>offset) continue;
		results[Npeaks].imageID = Nimage;
		results[Npeaks].pixelID = peak_y_tmp[i] * size_x + peak_x_tmp[i];
		results[Npeaks].peak_x = peak_x_tmp[i] + vars[0] - offset;
		results[Npeaks].peak_y = peak_y_tmp[i] + vars[1] - offset;
		results[Npeaks].intensity = vars[2] * ABscale;
		results[Npeaks].sigma_x = sigma_x;
		results[Npeaks].sigma_y = sigma_y;
		results[Npeaks].background = bg_fit * ABscale;
		results[Npeaks].max_x = peak_x_tmp[i];
		results[Npeaks].max_y = peak_y_tmp[i];
		results[Npeaks].Ncounts = sum_i;
		results[Npeaks].chi2 = 0.;
		for (j = 0; j < osize_sq; j++) results[Npeaks].chi2 += f[j] * f[j];
		results[Npeaks].chi2 /= (osize_sq - NVARS + 1);
		results[Npeaks++].lm_message = info;

		delete[] f; delete[] weights;
		delete[] icut; delete[] fjac; delete[] wa4;
	}


	//	delete[] peak_x_tmp; delete[] peak_y_tmp; 
	delete[] badpeaks; delete[] vars; 
	delete[] diag; delete[] qtf;
	delete[] wa1; delete[] wa2; delete[] wa3;
	delete[] ipvt; delete[] s_int; delete[] s_dbl;

	return 0;

}

/////////////////////////////////////////////////////////////////////////////////////////////

// calculate C and D so that C*coeft_s = D; see also paint_sigma

int paint_bg(int* image,			// one or several images
		int Nimages,			// number of images
		int size_x, int size_y,		// image size
		double* C, double* D)		// C*coeft_s = D: to be solved externally
{
  
  int i, j, x, y, k;
  int size1 = size_x*size_y;
  double s;

  // "mean" image
  int* I = new int[size1];
  for(j=0; j<size1; j++) I[j] = image[j];
  for(i=1; i<Nimages; i++)
    for(j=0; j<size1; j++) I[j] += image[i*size1+j];

  // F matrix: F[ji] = j-th point, i-th function
  double* F = new double[size1*6];
  j = 0;
  for (y=0; y<size_y; y++)
    for (x=0; x<size_x; x++) F[j++] = x*x;
  for (y=0; y<size_y; y++)
    for (x=0; x<size_x; x++) F[j++] = x*y;
  for (y=0; y<size_y; y++)
    for (x=0; x<size_x; x++) F[j++] = y*y;
  for (y=0; y<size_y; y++)
    for (x=0; x<size_x; x++) F[j++] = x;
  for (y=0; y<size_y; y++)
    for (x=0; x<size_x; x++) F[j++] = y;
  for (y=0; y<size_y; y++)
    for (x=0; x<size_x; x++) F[j++] = 1.;

  // C matrix
  for (i=0; i<6; i++)
    for (j=i; j<6; j++) {
      s = 0.;
      for (k=0; k<size1; k++) s += F[i*size1+k]*F[j*size1+k];
      C[i*6+j] = s;
      C[j*6+i] = s;
    }

  // D
  for (j=0; j<6; j++) {
    s = 0.;
    for (k=0; k<size1; k++) s += I[k]*F[j*size1+k]/(double)Nimages;
    D[j] = s;
  }
  
  delete[] I; delete[] F;
  return 0;
}


// calculate bg_surface and std from it

double paint_sigma(int* image,			// one or several images
		int Nimages,			// number of images
		int size_x, int size_y,		// size
		double* coefs)			// polynom coef-ts in order: x^2, xy, y^2, x, y, 1

{
  int x, y, i, j;
  double a = coefs[0], b = coefs[1], c = coefs[2], 
         d = coefs[3], e = coefs[4], f = coefs[5];

  // calculate bg_surface
  if (bg_surface!=NULL) delete[] bg_surface;
  bg_surface = new double[size_x*size_y];
  for (y=0; y<size_y; y++)
    for (x=0; x<size_x; x++)
      bg_surface[y*size_x+x] = f + x*(d+a*x) + y*(e+b*x+c*y);

  if (threshold_surface!=NULL) delete[] threshold_surface;
  threshold_surface = new int[size_x*size_y];
  last_threshold = -108;

  double s = 0., dif, sigma6;
  int size1 = size_x*size_y, n = size1*Nimages;

  // sigma: first run
  for (i=0; i<Nimages; i++)
    for (j=0; j<size1; j++) {
      dif = image[i*size1+j] - bg_surface[j];
      s += dif*dif;
    }

  // sigma: second run. Discard if > 6*sigma
  sigma6 = 6.*sqrt(s/(double)(n));
  s = 0.;
  for (i=0; i<Nimages; i++)
    for (j=0; j<size1; j++) {
      dif = image[i*size1+j] - bg_surface[j];
      if (dif > sigma6) n--;
      else s += dif*dif;
    }

  return sqrt(s/(double)(n));
}

int get_bg_surface(double* bg_out, int size_x, int size_y) // e.g. for testing
{
  for(int i=0; i<size_x*size_y; i++) bg_out[i] = bg_surface[i];
  return 0;
}


#endif

// cl /O2 /EHsc /MD /W3 /Zp1 /Fepaint.dll paint_analysis_bg.cpp *.c /link /dll /def:paint.def
// x64: cl /O2 /EHsc /MD /W3 /Fepaint.dll paint_analysis_bg.cpp *.c /link /dll /def:paint.def