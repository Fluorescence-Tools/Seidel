// compile everything:
// cl /O2 /EHsc /MD /TP /Fefits2x.dll fit23.cpp fit24.cpp fit25.cpp fit26.cpp i_lbfgs.cpp lbfgs.cpp modelf.cpp fsconv2.c twoIstar.cpp libs/ap.cpp /link /dll /export:fit23 /export:fit24 /export:fit25 /export:fit26
// mt /manifest fits2x.dll.manifest /outputresource:"fits2x.dll;#2"
// LabView data types

#include <stdlib.h>
#include <pshpack1.h>
#include <poppack.h>

typedef struct {
  int length;
  int data[1];
} LVI32Array;

typedef struct {
  int length;
  double data[1];
} LVDoubleArray;

typedef struct {
  LVI32Array** expdata;
  LVDoubleArray** irf;
  LVDoubleArray** bg;	// must be normalized outside!!!
  double dt;
  LVDoubleArray** corrections;
  LVDoubleArray** M;
} MParam;

typedef struct {
	LVDoubleArray** subimage;
	int osize;
	LVDoubleArray** M;
} MGParam;

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

void init_fact();
double loggammaf(double);
double wcm_p2s(int, double, double);
double Wcm_p2s(int*, double*, int);
double Wcm(int*, double*, int);
double W2DG(double*, double*, int, int);
double twoIstar_p2s(int*, double*, int);
double twoIstar(int*, double*, int);

// normalization
static double Sp, Ss, Bp, Bs;
static double Bexpected;	// expected <B> corresponding to the mean Bg signal
static int fixedrho = 0;
static int softbifl = 0;
static int p2s_twoIstar = 0;
static int firstcall = 1;
static double penalty = 0.;

static void normM_p2s (double* M, int Nchannels)
{
  int i;
  double s = 0.;

  for(i=0; i<Nchannels; i++) s += M[i]; 
  for(i=0; i<Nchannels; i++) M[i] *= Sp/s;

  s = 0.;
  for(i=Nchannels; i<2*Nchannels; i++) s += M[i];
  for(i=Nchannels; i<2*Nchannels; i++) M[i] *= Ss/s;
}

static void normM (double* M, int Nchannels)
{
  int i;
  double s = 0., Sexp = Sp + Ss;

  for(i=0; i<2*Nchannels; i++) s += M[i]; 
  for(i=0; i<2*Nchannels; i++) M[i] *= Sexp/s;

}
// if already normalized to s:
static void normM (double* M, double s, int Nchannels)
{
  int i;
  double Sexp = Sp + Ss;

  for(i=0; i<2*Nchannels; i++) M[i] *= Sexp/s;
}

// Gauss_analysis
double fit2DGaussian(double* vars, MGParam * p);
int model2DGaussian(double* vars, double* pM, int osize);
int modelTwo2DGaussian(double* vars, double* M, int osize);
int modelThree2DGaussian(double* vars, double* M, int osize);
double target2DGaussian(double* vars, void* pM);
int Gauss2D_analysis_Ani(double*, int, int, OptsCluster*, int, int&, int&, int&, int&, void**, int, int, int*, int*, int, MGParam*);