// gradient and jacobian: forward-difference, two-points, 4-points approximations

#include <math.h>

//////////////////////////////////// forward-difference /////////////////////////////////

int fjac1(void (*f)(double*, double*), double* x, int m, int n, double eps, double* fjac)
{

 int ij;
 double h, temp;
 double* fvec = new double[m];
 double* wa = new double[m];

 f(x,fvec);

 ij = 0;
 for(int j=0; j<n; j++ ) {
   temp = x[j];
   h = eps * abs(temp);
   if(h == 0.) h = eps;
   x[j] = temp + h;
   f(x,wa);
   x[j] = temp;

   for(int i=0; i<m; i++ ) {
     fjac[ij] = (wa[i] - fvec[i])/h;
     ij += 1;	/* fjac[i+m*j] */
   }
 }

 delete[] fvec; delete[] wa; 
 return 0;
}

int fgrad1(void (*f)(double*, double&), double* x, int n, double eps, double* fgrad)
{

 double h, temp;
 double fval;
 double w;

 f(x,fval);

 for(int j=0; j<n; j++ ) {
   temp = x[j];
   h = eps * abs(temp);
   if(h == 0.) h = eps;
   x[j] = temp + h;
   f(x,w);
   x[j] = temp;

   fgrad[j] = (w - fval)/h;
 }

 return 0;
}


///////////////////////////////////////// 2 points //////////////////////////////////////

int fjac2(void (*f)(double*, double*), double* x, int m, int n, double eps, double* fjac)
{

 int ij;
 double h, temp;
 double* wa1 = new double[m];
 double* wa2 = new double[m];

 ij = 0;
 for(int j=0; j<n; j++ ) {
   temp = x[j];
   h = eps * abs(temp);
   if(h == 0.) h = eps;
   x[j] = temp + h;
   f(x,wa1);
   x[j] = temp - h;
   f(x,wa2);
   x[j] = temp;

   for(int i=0; i<m; i++ ) {
     fjac[ij] = 0.5*(wa1[i] - wa2[i])/h;
     ij += 1;	/* fjac[i+m*j] */
   }
 }

 delete[] wa1; delete[] wa2; 
 return 0;
}

int fgrad2(void (*f)(double*, double&), double* x, int n, double eps, double* fgrad)
{

 double h, temp;
 double w1, w2;

 for(int j=0; j<n; j++ ) {
   temp = x[j];
   h = eps * abs(temp);
   if(h == 0.) h = eps;
   x[j] = temp + h;
   f(x,w1);
   x[j] = temp - h;
   f(x,w2);
   x[j] = temp;

   fgrad[j] = 0.5*(w1 - w2)/h;
 }

 return 0;
}

///////////////////////////////////////// 4 points //////////////////////////////////////

int fjac4(void (*f)(double*, double*), double* x, int m, int n, double eps, double* fjac)
{

 const double c1 = 2./3.;
 const double c2 = 1./12.;

 int ij;
 double h, temp;
 double* wa1 = new double[m];
 double* wa2 = new double[m];

 ij = 0;
 for(int j=0; j<n; j++ ) {
   temp = x[j];
   h = eps * abs(temp);
   if(h == 0.) h = eps;
   x[j] = temp + h;
   f(x,wa1);
   x[j] = temp - h;
   f(x,wa2);

   for(int i=0; i<m; i++ ) {
     fjac[ij] = c1*(wa1[i] - wa2[i])/h;
     ij += 1;	/* fjac[i+m*j] */
   }

   ij -= m;
   x[j] = temp + 2.*h;
   f(x,wa1);
   x[j] = temp - 2.*h;
   f(x,wa2);
   x[j] = temp;

   for(int i=0; i<m; i++ ) {
     fjac[ij] += c2*(wa2[i] - wa1[i])/h;
     ij += 1;	/* fjac[i+m*j] */
   }

 }

 delete[] wa1; delete[] wa2; 
 return 0;
}

int fgrad4(void (*f)(double*, double&), double* x, int n, double eps, double* fgrad)
{

 const double c1 = 2./3.;
 const double c2 = 1./12.;

 double h, temp;
 double w1, w2;

 for(int j=0; j<n; j++ ) {
   temp = x[j];
   h = eps * abs(temp);
   if(h == 0.) h = eps;
   x[j] = temp + h;
   f(x,w1);
   x[j] = temp - h;
   f(x,w2);

   fgrad[j] = c1*(w1 - w2)/h;

   x[j] = temp + 2.*h;
   f(x,w1);
   x[j] = temp - 2.*h;
   f(x,w2);
   x[j] = temp;

   fgrad[j] += c2*(w2 - w1)/h;
 }

 return 0;
}