// minimize f(x,p) using BFGS algorithm

#include <math.h>
#include "lbfgs.h"

// pointer to the target function
typedef double(*TargetFP)(double*, void*);

#ifdef WIN32
#pragma pack(1)
#endif
class bfgs
{

 public:
  
  bfgs(TargetFP fun) { setdefaults(); N = 0; f = fun; }
  bfgs(TargetFP fun, int n) { setdefaults(); setN(n); f = fun; }
  ~bfgs() {
    if (N>0) { delete[] xd; delete[] fixed; }
  }
  
  // change dimension
  void setN(int n) {
    if (N>0) { delete[] xd; delete[] fixed; }
    N = n;
    xd = new double[N];
    fixed = new int[N];
    for(int i=0; i<N; i++) fixed[i] = 0;
  }

  // set epsilon explicitly
  void seteps(double e) {
    eps = e;
    sqrt_eps = sqrt(e);
  }
  
  // estimate epsilon
  void seteps() {
    double x1 = 1.0, x2 = 1.0, e = 1.e-18, estep = 1.001;
    for (;;) {
      x2 = x1 + e;
      if (x2 > x1) break;
      else e *= estep;
    }
    seteps(3.*e);
  }

  // fix or unfix a paramter
  void fix(int n) {
    if (n>=N) return;
    else fixed[n] = 1;
  }
  void free(int n) {
    if (n>=N) return;
    else fixed[n] = 0;
  }

  // max number of iterations
  int maxiter;
  
  // minimize target f
  int minimize(double* x, void* p)
  {

    if (N==0) return -1;

    int info = 0, bfgsM, Nfree = 0, j = 1;
    for(int i=0; i<N; i++) {
      xd[i] = x[i];
      Nfree += (!fixed[i]);
    }
    Nfree < 7 ? bfgsM = Nfree : bfgsM = 7;

    // create "real_1d_array" X for lbfgsminimize
    ap::real_1d_array X;
    X.setbounds(1, Nfree);
    for(int i=0; i<N; i++) 
      if (!fixed[i]) X(j++) = x[i];

    pcopy = p;
    lbfgsminimize(Nfree, bfgsM, X, 100.*sqrt_eps, 10000.*eps, 100.*sqrt_eps, maxiter, info);

    // copy results back to x
    j = 1;
    for(int i=0; i<N; i++) 
      if (!fixed[i]) x[i] = X(j++);

    return info;
  }

 private:

  int N;
  double eps;
  double sqrt_eps;
  TargetFP f;
  void* pcopy;
  double* xd;
  int* fixed;
  
  void setdefaults()
  {
     N = 0;
     // seteps();
	 seteps(3.3e-16);
     maxiter = 100;
  }

  // subroutines from lbfgs
  void lbfgsminimize(const int&, const int&, ap::real_1d_array&, const double&, const double&,
     const double&, const int&, int&);
  void lbfgslincomb(const int&, const double&, const ap::real_1d_array&, 
     int, ap::real_1d_array&, int);
  double lbfgsdotproduct(const int&, const ap::real_1d_array&, int, 
     const ap::real_1d_array&, int);
  void lbfgsmcsrch(const int&, ap::real_1d_array&, double&, ap::real_1d_array&,
     const ap::real_1d_array&, int, double&, const double&, const double&, const int&,
     int&, int&, ap::real_1d_array&, const double&, const double&, const double&);
  void lbfgsmcstep(double&, double&, double&, double&, double&, double&, double&,
     const double&, const double&, bool&, const double&, const double&, int&);
  void lbfgsnewiteration(const ap::real_1d_array&, double, const ap::real_1d_array&);

  // f gradient: 1, 2, and 4-points approximations
  void fgrad1(ap::real_1d_array&, double&, ap::real_1d_array&);
  void fgrad2(ap::real_1d_array&, double&, ap::real_1d_array&);
  void fgrad4(ap::real_1d_array&, double&, ap::real_1d_array&);
  void funcgrad(ap::real_1d_array& x, double& fval, ap::real_1d_array& g) { fgrad2(x, fval, g); }

};