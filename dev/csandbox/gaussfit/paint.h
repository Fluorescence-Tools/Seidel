// header for paint_analysis, paint_bg, ...

#include <stdlib.h>

// fitting options
typedef struct {
	int threshold; 			// threshold for peak search
	char relative_thr;		// threshold is relative to bg_surface
	int offset;			// offset (max +/- deltax, deltay) for fitting
	int maxNPeaks;
	char fit2DGauss;		// fit or just return cm's
	double sigma0;			// initial estimation for sigma
	double dsigma_max;		// threshold for |sigma_x-sigma_y|/(sigma_x+sigma-y)
	char must_converge;		// discard peaks for which LM has not converged
	int max_threshold; 			// max_threshold for peak search
		} OptsCluster;

// fit results
typedef struct {
	unsigned int imageID;
	unsigned int pixelID;
	double peak_x;
	double peak_y;
	double intensity;
	double chi2;
	int lm_message;
	double std;
	double sigma_x;
	double sigma_y;
	double background;
	double max_x;
	double max_y;	
	double Ncounts;
	} ResultsCluster;

// LM optimisation procedure

extern "C" int lmdif(double*, double*, double*, int, int, int,
          double*, double*, int*, double*, double*, double*, double*, double*,
          int*, double*, int*, int*);


// paint_analysis
void wr_2DGauss(double*, double*, int, double*, double*, double);
void wr_2DGauss_grad(double*, double*, int, double*, double);
int paint_analysis(int*, int, int, OptsCluster*, int, int&, void**);
int Gauss2D_analysis(int*, int, int, OptsCluster*, int, int&, void**, int, int*, int*);
