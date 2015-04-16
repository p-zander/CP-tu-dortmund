#include <functional>
#include <iostream>
#include <cmath>

using namespace std;

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

void findNforEps(double bottom, double top, unsigned int N, double eps,
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
  cout << "N: " << N << "	" << result1 << endl;
}

int main() {
  findNforEps(1, 100, 100, 1e-4, [](double x) { return exp(-x) / x; }, simpson);
  findNforEps(1, 100, 100, 1e-4, [](double x) { return exp(-x) / x; }, trapez);
  findNforEps(1, 100, 100, 1e-4, [](double x) { return exp(-x) / x; }, mp);

  // set lower limit to 1e-20 so division doesn't throw an error.
  findNforEps(1e-20, 1, 100, 1e-4, [](double x) { return x * sin(1 / x); },
              simpson);
  findNforEps(1e-20, 1, 100, 1e-4, [](double x) { return x * sin(1 / x); },
              trapez);
  findNforEps(1e-20, 1, 100, 1e-4, [](double x) { return x * sin(1 / x); }, mp);

  return 0;
}
