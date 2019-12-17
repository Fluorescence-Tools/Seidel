// weighted residuals with and without scaling for sfit/gfit program
// Version 2007-11-17

#include <math.h>
#include "fsconv.h"
#include "mex.h"

/* generic weighted residuals */
void wres(double *wr, double *fit, double *decay, double *w_sq, 
          int start, int stop)
{
  int i;
  double w;

  for (i=start; i<=stop; i++) {
    w = sqrt(w_sq[i]);
    if (w==0) w = 1.;
    wr[i-start] = (decay[i]-fit[i])/w;
  }
}

/* sconv via MatLab sconv_fft */
void sconv_matlab(double *fit, double *p, double *lamp, int start, int stop)
{
 int i;
 mxArray *sconv_fft_in[4], *sconv_fft_out[1];
 char* sconv_fft_name = "sconv_fft";
 double* mxpointer;

 // model function
 sconv_fft_in[0] = mxCreateDoubleMatrix(stop+1, 1, mxREAL); 
 mxpointer = mxGetPr(sconv_fft_in[0]);
 for (i=0; i<=stop; i++) mxpointer[i] = p[i];

 // lamp
 sconv_fft_in[1] = mxCreateDoubleMatrix(stop+1, 1, mxREAL);
 mxpointer = mxGetPr(sconv_fft_in[1]);
 for (i=0; i<=stop; i++) mxpointer[i] = lamp[i];

 // fit range
 sconv_fft_in[2] = mxCreateDoubleMatrix(1, 2, mxREAL);
 mxpointer = mxGetPr(sconv_fft_in[2]);
 mxpointer[0] = 1.; mxpointer[1] = (double)(stop+1);

 // DELTA_T
 sconv_fft_in[3] = mxCreateDoubleMatrix(1, 1, mxREAL);
 mxpointer = mxGetPr(sconv_fft_in[3]);
 mxpointer[0] = DELTA_T;

 mexCallMATLAB(1, sconv_fft_out, 4, sconv_fft_in, sconv_fft_name);

 // return
 mxpointer = mxGetPr(sconv_fft_out[0]);
 for (i=start; i<=stop; i++) fit[i] = mxpointer[i];

 mxDestroyArray(sconv_fft_in[0]);
 mxDestroyArray(sconv_fft_in[1]);
 mxDestroyArray(sconv_fft_in[2]);
 mxDestroyArray(sconv_fft_in[3]);
 mxDestroyArray(sconv_fft_out[0]);

}


void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[])
  
{ 
 double *xs, *xd, *lamp, *corr, *decay, *vv, *vh, *range;
 double *wr;
 double lampsh[N_CHANNELS];
 double ts, bg, bg_vv, bg_vh, scatter, period, g, l1, l2, scale=0., s_orig;
 double sumsum, sumdif, difsum, difdif;
 int vm, len_xs, len_xd, fitstart, fitstop, i, n_points;
 int autoscale, fit_vv_vh;
 double fit_s[N_CHANNELS], fit_d[N_CHANNELS], curve[N_CHANNELS], w2[N_CHANNELS];

/************************** Check input arguments **************************/

 /* Arguments allowed:
        1) x = [a1 tau1 a2 tau2 ...],
        corrections = [ts bg(_vv) {bg_vh} scatter autoscale period {g fit_vv_vh l1 l2}],
        lamp, decay, range = [fitstart fitstop];
        2) xsum, xdif, ..., vv, vh, ...
        or: p(t) arrays instead of x*/

 if (nrhs==5) vm=1;
 else if (nrhs==7) vm=0;
 else if (nrhs==0) {	// return DELTA_T and N_CHANNELS
   plhs[0] = mxCreateDoubleMatrix(1, 2, mxREAL);
   wr = mxGetPr(plhs[0]);
   wr[0] = N_CHANNELS;
   wr[1] = DELTA_T;
   return;
 }
 else {
   mexPrintf("Current settings: # channels = %d, dt = %f\n", N_CHANNELS, DELTA_T);
   mexErrMsgTxt("Error: 5 (lifetime experiment) or 7 (anisotropy experiment)"
 	" argumets allowed.");
 }

/************************ Input and output arguments ***********************/

 xs = mxGetPr(prhs[0]);
 len_xs = mxGetN(prhs[0])*mxGetM(prhs[0]);

 if (vm) {
  corr = mxGetPr(prhs[1]);
  lamp = mxGetPr(prhs[2]);
  n_points = mxGetN(prhs[2])*mxGetM(prhs[2]);
  decay = mxGetPr(prhs[3]);
  range = mxGetPr(prhs[4]);
 }

 else { /* not vm */
  xd = mxGetPr(prhs[1]);
  len_xd = mxGetN(prhs[1])*mxGetM(prhs[1]);
  corr = mxGetPr(prhs[2]);
  lamp = mxGetPr(prhs[3]);
  n_points = mxGetN(prhs[3])*mxGetM(prhs[3]);
  vv = mxGetPr(prhs[4]);
  vh = mxGetPr(prhs[5]);
  range = mxGetPr(prhs[6]);
 }

 ts = corr[0]/DELTA_T;
 if (vm) {
   bg = corr[1];
   scatter = corr[2];
   autoscale = (int)corr[3];
   period = corr[4];
   }
 else { // !vm
   bg_vv = corr[1];
   bg_vh = corr[2];
   scatter = corr[3];
   autoscale = (int)corr[4];
   period = corr[5];
   g = corr[6];
   fit_vv_vh = (int)corr[7];
   l1 = corr[8];
   l2 = corr[9];
 }
 fitstart = (int)range[0]-1; fitstop = (int)range[1]-1;

 /* Create matrix and pointers for the return arguments */ 
 plhs[0] = mxCreateDoubleMatrix(fitstop-fitstart+1, 2-vm, mxREAL); 

 shift_lamp(lampsh,lamp,ts,n_points);
 wr = mxGetPr(plhs[0]);

/*************************** Weigths and residuals *************************/

 /* convolution */
 if (len_xs<fitstop+1) {
   fconv_per(fit_s, xs, lampsh, len_xs/2, fitstart, fitstop, n_points, period);
   if (!vm) fconv_per(fit_d, xd, lampsh, len_xd/2, fitstart, fitstop, n_points, period);
 }
 else {
   sconv_matlab(fit_s, xs, lampsh, fitstart, fitstop);
   if (!vm) sconv_matlab(fit_d, xd, lampsh, fitstart, fitstop);
 }

 /* scatter */
 for (i=fitstart; i<=fitstop; i++) fit_s[i] += scatter*lamp[i];
 if (!vm)  for (i=fitstart; i<=fitstop; i++) fit_d[i] += scatter*lamp[i];

 /* l1, l2 */
 if (!vm) {
   sumsum = l1/3.*(-1.+1./g) + l2/3.*(2.*g-2.);
   sumdif = l1/3.*(-1.+1./g) + l2/3.*(-g+1.);
   difsum = l1/3.*(-2.-1./g) + l2/3.*(4.*g+2.);
   difdif = l1/3.*(-2.-1./g) + l2/3.*(-2.*g-1.);
   for (i=fitstart; i<=fitstop; i++) {
     s_orig = fit_s[i];
     fit_s[i] += sumsum*s_orig + difsum*fit_d[i];
     fit_d[i] += sumdif*s_orig + difdif*fit_d[i];
   }
 }     

 if (vm) {

   if (autoscale) rescale_w_bg(fit_s, decay, decay, bg, &scale, fitstart, fitstop);
   for (i=fitstart; i<=fitstop; i++) fit_s[i] += bg;
   wres(wr, fit_s, decay, decay, fitstart, fitstop);
 }

 else { /* not vm */  

  for (i=fitstart; i<=fitstop; i++) {
   w2[i] = vv[i] + 4.*g*g*vh[i];
   curve[i] = vv[i] + 2.*g*vh[i];
  }
  if (autoscale) {
    rescale_w_bg(fit_s, curve, w2, bg_vv + 2.*g*bg_vh, &scale, fitstart, fitstop);
    rescale_w_bg(fit_d, curve, w2, bg_vv - g*bg_vh, &scale, fitstart, fitstop);
  }
  for (i=fitstart; i<=fitstop; i++) {
    fit_s[i] += bg_vv + 2.*g*bg_vh;
    fit_d[i] += bg_vv - g*bg_vh;
  }

  if (fit_vv_vh) {
    for (i=fitstart; i<=fitstop; i++)  curve[i] = (fit_s[i] + 2.*fit_d[i])/3.;
    wres(wr, curve, vv, vv, fitstart, fitstop);
    for (i=fitstart; i<=fitstop; i++)  curve[i] = (fit_s[i] - fit_d[i])/3./g;
    wres(wr+(fitstop-fitstart+1), curve, vh, vh, fitstart, fitstop);
  }
  else { /* fit sum and dif */
    wres(wr, fit_s, curve, w2, fitstart, fitstop);
    for (i=fitstart; i<=fitstop; i++) {
      w2[i] = vv[i] + g*g*vh[i];
      curve[i] = vv[i] - g*vh[i];
    }
    wres(wr+(fitstop-fitstart+1), fit_d, curve, w2, fitstart, fitstop);
  }

 }

 return;
}
