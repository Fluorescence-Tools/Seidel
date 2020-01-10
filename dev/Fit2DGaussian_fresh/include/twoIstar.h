#ifndef TWOISTAR_
#define TWOISTAR_
/*! return weighted MLE based on data C, model M of length Nchannels
	bool wi_nowi is deprecated
*/
double W2DG(double* data, double* model, int osize);
double twoIstar_G(double* C, double* M, int osize);

#endif