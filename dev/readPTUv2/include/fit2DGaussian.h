#ifndef FIT2DGAUSSIAN
#define FIT2DGAUSSIAN


#include <stdlib.h>
#include <pshpack1.h>
#include <poppack.h>

//Instead, use GaussDataType

typedef struct {
	double* data; //contains the data to be fitted
	double* model;//initialise empty array that will contain the model
	int xlen;
	int ylen; // for 1D data, ylen is unused
} GaussDataType;

// fitting options
typedef struct {
	char elliptical_circular;	// circular ?
	double background; 			// const BG input
	char free_fixed;			// with free or fixed BG ?
	int maxNPeaks;
	char fit2DGauss;			// fit or just return cm's
	char must_converge;			// discard peaks for which LM has not converged
} OptsCluster;

// fit results
// legacy from Suren
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

#endif
