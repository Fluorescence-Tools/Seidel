#ifndef EIGEN_WRAP
#define EIGEN_WRAP

//get import needed in function declarations
#include <Eigen/Core>
using namespace Eigen;


/*! inverse a numpy matrix
*/
Eigen::MatrixXd inv(Eigen::MatrixXd xs);

/*! Substract one integer from another
	\param i an integer
	\param j an integer to subtract from \p i
*/
int subtract(int i, int j);

#endif