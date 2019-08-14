// 2I* parameter for Sp + 2Ss

#include<math.h>
#include <iostream>
using namespace std;
#include <time.h>

// factorial
static double fact[150];
void init_fact()
{
  double f = 1.;
  fact[0] = 1.;
  for(int i = 1; i<150; i++) {
    f *= (double)i;
    fact[i] = f;
  }
}

// log-likelihood w(C|m) for Cp + 2Cs

double wcm(int C, double m)
{
  return C*log(m);
}

// overall log-likelihood w(C,M)
double twoIstar_p2s(int* C, double* M, int Ndata)
{
  double W = 0., W0 = 0.;
  int nempty = 0;
  for (int i=0; i<Ndata; i++)
    if (C[i]>0) {
      W += wcm(C[i], M[i]);
      W0 += wcm(C[i], (double)C[i]);
    }
    else {W += 1.; W0 += 1.;} // nempty++;
  return -2.*(W-W0)/(double)Ndata;
}

int main()
{
  srand((unsigned)time( NULL ));
  init_fact();
  const int Ndata = 1024;
  int nt, ntbest, repeats = 100000;
  int* d = new int[Ndata];
  double* M = new double[Ndata];
  double* thist = new double[41];
  double s = 0., sm = 0., smean = 0., twoIbest, twoI, taubest;

  for (int i=0; i<41; i++) thist[i] = 0.;

for (int j = 0; j<repeats; j++) {
  for (int i=0; i<Ndata; i++) d[i] = 0.;
  s = 0.;
  double t, dt = 52./Ndata;
  for (int i=0; i<200; i++) {
    t = -4.*log(1.-rand()/(double)RAND_MAX);
    if (t/dt < (double)Ndata) { d[(int)(t/dt)]++; s++; }
  }

  smean += s/(double)repeats;
  twoIbest = 100.;
  nt = 0;
/*  for (t=2; t<6; t += 0.1) {
    sm = 0.;
    for (int i=0; i<Ndata; i++) { M[i] = exp(-i*dt/t); sm += M[i]; }
    for (int i=0; i<Ndata; i++) { M[i] *= s/sm; }
    twoI = twoIstar_p2s(d, M, Ndata);
    if (twoI < twoIbest) {twoIbest = twoI; taubest = t; ntbest=nt;}

    nt++;
  }
   thist[ntbest]++; */
   sm = 0.;
   for (int i=0; i<Ndata; i++) { sm += (i*dt)*d[i]; }

   ntbest = (int)(10.*(sm/s -2.));
   if (ntbest > 40) ntbest = 40;
   thist[ntbest]++;
}
  nt = 0;
  for (double t=2; t<6; t += 0.1) {
   cout << thist[nt] << endl;
   nt++;
  }
  // chi2
  double chi2 = 0.;
  for (int i=0; i<Ndata; i++) chi2 += (d[i] - M[i])*(d[i] - M[i])/(M[i]);
  cout << "chi2 = " << chi2 / (double)Ndata << endl;
  cout << "<s> = " << smean << endl;

  return 1;
}
