// model function for fits 23, 24, ...
// Version 2008-04-24

#include <math.h>
#include "fsconv.h"

/////////////////////////////////////////// fit23 /////////////////////////////////////////

int modelf23(double* param, 			// here: [tau gamma r0 rho]
		double* irf,
		double* bg,
		int Nchannels,
		double dt,			// time per channel
		double* corrections,		// [period g l1 l2]
		double* mfunction)		// out: model function in Jordi-girl format
  
{ 
 double x[4];
 double tau, gamma, r0, rho, 
  period, g, l1, l2, const_bg,
  taurho, sum_m = 0., sum_s = 0.,
  f1, f3;
 int i, conv_stop;

/************************ Input arguments ***********************/

 tau = param[0]*DELTA_T/dt; gamma = param[1];
 r0 = param[2]; rho = param[3]*DELTA_T/dt;

 period = corrections[0]*DELTA_T/dt; g = corrections[1];
 l1 = corrections[2]; l2 = corrections[3];
 conv_stop = (int)corrections[4];
 const_bg = 0.01 * corrections[5];

/************************* Model function ***********************/

 taurho = 1./(1./tau + 1./rho);

 /// vv + vh

 x[0] = 1.; x[1] = tau;
 x[2] = r0 * (2.-3.*l1); x[3] = taurho;
 fconv_per_cs(mfunction, x, irf, 2, Nchannels-1, Nchannels, period, conv_stop);
 
 x[0] = 1./g; x[2] = 1./g * r0 * (-1.+3.*l2);
 fconv_per_cs(mfunction + Nchannels, x, irf + Nchannels, 2, Nchannels-1, Nchannels, period, conv_stop);
 
 /// add scatter

 for (i=0; i<2*Nchannels; i++) sum_m += mfunction[i]; 
 f1 = (1.-gamma-const_bg)/sum_m; f3 = 0.5*const_bg/(double)Nchannels;
 for (i=0; i<2*Nchannels; i++)
   mfunction[i] = mfunction[i]*f1 + bg[i]*gamma + f3;

 return 0;
}

/////////////////////////////////////////// fit24 /////////////////////////////////////////

int modelf24(double* param, 			// here: [tau1 gamma tau2 A2 offset]
		double* irf,
		double* bg,
		int Nchannels,
		double dt,			// time per channel
		double* corrections,		// [period g l1 l2]
		double* mfunction)		// out: model function in Jordi-girl format
  
{ 
 double x[4];
 double tau1, gamma, tau2, A2, offset, 
  period,
  sum_m = 0., sum_s = 0.;
 int i, conv_stop;

/************************ Input arguments ***********************/

 tau1 = param[0]*DELTA_T/dt; gamma = param[1];
 tau2 = param[2]*DELTA_T/dt; A2 = param[3];
 offset = param[4]/(double)Nchannels;

 period = corrections[0]*DELTA_T/dt;
 conv_stop = (int)corrections[4];

/************************* Model function ***********************/

 /// vv

 x[0] = 1.-A2; x[1] = tau1;
 x[2] = A2; x[3] = tau2;
 fconv_per_cs(mfunction, x, irf, 2, Nchannels-1, Nchannels, period, conv_stop);

 /// vh

 fconv_per_cs(mfunction + Nchannels, x, irf + Nchannels, 2, Nchannels-1, Nchannels, period, conv_stop);

 /// add scatter and background

 for (i=0; i<2*Nchannels; i++) { sum_m += mfunction[i]; sum_s += bg[i]; }
 for (i=0; i<2*Nchannels; i++)
   mfunction[i] = mfunction[i] * (1.-gamma)/sum_m + bg[i]*gamma/sum_s + offset;

 return 0;

}
