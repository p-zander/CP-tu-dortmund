
#include "iostream"
#include "functional"
#include "cmath"

using std::cout; using std::endl; using std::function; using std::abs;

double trapez(double bottom, double top, unsigned int N,
              function<double(double)> const& f) {
  double h = (top - bottom) / N;
  double sum = 0;

  for (unsigned int a = 1; a < N; ++a) {
    sum += f(bottom + a * h);
  }
  sum *= h;
  sum += h / 2 * (f(bottom) + f(top));
  return sum;
}

double simpson(double bottom, double top, unsigned int N,
               function<double(double)> const& f) {
  double h = (top - bottom) / N;
  double sum = 0;
  N = (N % 2) ? N = 1 : N;
  for (unsigned int a = 1; a < N; ++a) {
    sum += (a % 2) ? 4 * f(bottom + a * h) : 2 * f(bottom + a * h);
  }
  sum += f(bottom) + f(top);
  sum *= h / 3;
  return sum;
}

double mp(double bottom, double top, unsigned int N,
          function<double(double)> const& f) {
  double h = (top - bottom) / N;
  double sum = 0;
  for (unsigned int a = 1; a < N; ++a) {
    sum += f(bottom + (a + 0.5) * h);
  }
  sum *= h;
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
	//cout << "N: " << N << "	" << result1 << endl;

  return result1;
}

int main(int argc, char const *argv[]) {
	double simpsonResult = findNforEps(-1, 1, 100, 1e-7, [](double x) { if(x==0) return 1.; return ( exp(x) - 1 )/ x; }, simpson);
	double  trapezResult = findNforEps(-1, 1, 100, 1e-7, [](double x) { if(x==0) return 1.; return ( exp(x) - 1 )/ x; }, trapez);
	double      mpResult = findNforEps(-1, 1, 100, 1e-7, [](double x) { if(x==0) return 1.; return ( exp(x) - 1 )/ x; }, mp);

	cout.precision(10);

	cout << "I_1:" << endl << "Simpson: " << simpsonResult << " , Trapez: " << trapezResult << " , MP: " << mpResult << endl;

	double mpResult_1 = findNforEps(-1, 0, 100, 1e-7, [](double x) { if(x==0) return 1.; return ( exp(x) - 1 )/ x; }, mp);
	double mpResult_2 = findNforEps( 0, 1, 100, 1e-7, [](double x) { if(x==0) return 1.; return ( exp(x) - 1 )/ x; }, mp);

	cout << "MP left+right: " << mpResult_1+mpResult_2 << endl;

	double I_2 = 2*findNforEps( 0, 1e7, 1e9, 1e-7, [](double x){ return exp(-x*x); }, mp);

	cout << "I_2 = " << I_2 << endl;

	return 0;
}