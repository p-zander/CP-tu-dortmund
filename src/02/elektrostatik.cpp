#include "cmath"
#include "iostream"
#include "functional"

using std::cout;
using std::endl;
using std::function;

double simpson(double bottom, double top, unsigned int N, function<double(double)> const& f) {
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

int main(int argc, char const* argv[]) {
    double a = 1;
    int steps = 100;

    for (int n = 1; n <= 80; n++) {
        if (n==10) continue;
        double X = n * 0.1 * a;
        double result = simpson(-a, a, steps, [X, a, steps](double x) {
                            return simpson(-a, a, steps, [X, a, steps, x](double y) {
                                return simpson(-a, a, steps, [X, a, steps, x, y](double z) {
                                    if ( X == x && y == 0 && z == 0 ) return 0.;
                                    return a / sqrt((X - x) * (X - x) + y * y + z * z);
                                });
                            });
                        });
        cout << X << "\t" << result << endl;
    }

    // cout << "Result at X = " << X << " is " << result << endl;
    return 0;
}
