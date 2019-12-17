// minimize f(x,p) using BFGS algorithm

#include "i_lbfgs.h"

////////////////////////////////////// gradient ////////////////////////////////////

void bfgs::fgrad1(ap::real_1d_array& x, double& fval, ap::real_1d_array& fgrad)
{

 double h, temp, w;
 int j = 1;

 for(int i=0; i<N; i++) 
   if (!fixed[i]) xd[i] = x(j++);

 fval = f(xd, pcopy);

 j = 1;
 for(int i=0; i<N; i++) {
   if (fixed[i]) continue;
   temp = xd[i];
   h = sqrt_eps * fabs(temp);
   if(h == 0.) h = sqrt_eps;
   xd[i] = temp + h;
   w = f(xd, pcopy);
   xd[i] = temp;

   fgrad(j++) = (w - fval)/h;
 }

}

///////////////////////////////////////// 2 points //////////////////////////////////////

void bfgs::fgrad2(ap::real_1d_array& x, double& fval, ap::real_1d_array& fgrad)
{

 double h, temp, w1, w2;
 int j = 1;

 for(int i=0; i<N; i++) 
   if (!fixed[i]) xd[i] = x(j++);

 fval = f(xd, pcopy);

 j = 1;
 for(int i=0; i<N; i++) {
   if (fixed[i]) continue;
   temp = xd[i];
   h = sqrt_eps * fabs(temp);
   if(h == 0.) h = sqrt_eps;
   xd[i] = temp + h;
   w1 = f(xd, pcopy);
   xd[i] = temp - h;
   w2 = f(xd, pcopy);
   xd[i] = temp;

   fgrad(j++) = 0.5*(w1 - w2)/h;
 }

}

///////////////////////////////////////// 4 points //////////////////////////////////////

void bfgs::fgrad4(ap::real_1d_array& x, double& fval, ap::real_1d_array& fgrad)
{

 const double c1 = 2./3.;
 const double c2 = 1./12.;

 double h, temp, w1, w2;
 int j = 1;

 for(int i=0; i<N; i++) 
   if (!fixed[i]) xd[i] = x(j++);

 fval = f(xd, pcopy);

 j = 1;
 for(int i=0; i<N; i++) {
   if (fixed[i]) continue;
   temp = xd[i];
   h = sqrt_eps * fabs(temp);
   if(h == 0.) h = sqrt_eps;
   xd[i] = temp + h;
   w1 = f(xd, pcopy);
   xd[i] = temp - h;
   w2 = f(xd, pcopy);

   fgrad(j) = c1*(w1 - w2)/h;

   xd[i] = temp + 2.*h;
   w1 = f(xd, pcopy);
   xd[i] = temp - 2.*h;
   w2 = f(xd, pcopy);
   xd[i] = temp;

   fgrad(j++) += c2*(w2 - w1)/h;
 }
}