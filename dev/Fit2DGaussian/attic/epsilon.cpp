#include <iostream>
using namespace std;

int main()
{
  double x1 = 1.0, x2 = 1.0;
  double e = 1.e-18;
  const double step = 1.001;
  for (;;) {
    x2 = x1 + e;
    if (x2 > x1) break;
    else e *= step;
  }
  cout << "epsilon = " << e << endl;
}