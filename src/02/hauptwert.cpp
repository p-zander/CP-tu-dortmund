#include "iostream"
#include "functional"
#include "cmath"

using std::cout; using std::endl; using std::function; using std::abs;

double simpson(double bottom, double top, unsigned int N,
               function<double(double)> const& f) {
  double h = (top - bottom) / N;
  double sum = 0;
  N = (N % 2) ? 1 : N;
  for (unsigned int a = 1; a < N; ++a) {
    sum += (a % 2) ? 4 * f(bottom + a * h) : 2 * f(bottom + a * h);
  }
  sum += f(bottom) + f(top);
  sum *= h / 3;
  return sum;
}

double findNforEps(double bottom, double top, unsigned int N, double eps,
                    function<double(double)> const& f,
                    function<double(double, double, unsigned int,
                    function<double(double)> const& f)> const& g) {
  double error = 1;
  double result1 = g(bottom, top, N, f);
  double result2 = 0;

  do {
    result2 = result1;
    N *= 2;
    result1 = g(bottom, top, N, f);
    error = abs(result2 - result1) / result2;
  } while (error > eps);

  return result1;
}

int main() {
    double I_1 =   findNforEps(-1,   1, 100, 1e-9, [](double x){ if(x==0) return 1.; return ( exp(x) - 1 )/ x; }, simpson);
    double I_2 = 2*findNforEps( 0, 100, 100, 1e-9, [](double x){ return exp(-x*x); }, simpson);

    cout.precision(10);
    cout << "I_1 = " << I_1 << " (Abweichung: " << I_1-2.11450175 << ")" <<endl;
    cout << "I_2 = " << I_2 << " (Abweichung: " << I_2-sqrt(M_PI) << ")" <<endl;

    return 0;
}
