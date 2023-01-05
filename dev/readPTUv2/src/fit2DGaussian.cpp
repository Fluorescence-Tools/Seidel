// Fit 2DGauss
// Version 2019-04-16
//edited extensively by Nicolaas van der Voort
//AG Seidel, Düsseldorf

#include <iostream>
#include <fstream>
#include <math.h>
//#include <set>
#include <Eigen/Core>

using namespace std;

const int NVARS = 6;
const int NPEAKS_FACTOR = 10;

typedef Eigen::Matrix<double, 5, 1> Vector5d;
typedef Eigen::Matrix<double, 6, 1> Vector6d;

// background and threshold surfaces
static double* bg_surface = NULL;
static int* threshold_surface = NULL;
static int last_threshold = -108;

int osize;//object size
int osize_sq;//object area
int* icut;

