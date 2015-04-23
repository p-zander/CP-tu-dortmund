#include "cmath"
#include "fstream"
#include "functional"

using std::endl;
using std::function;
using std::ofstream;

double simpson(double bottom, double top, unsigned int N, function<double(double)> const& f) {
    double h = (top - bottom) / N;
    double sum = 0;
    N = (N % 2) ? 1 : N;
    for (unsigned int a = 1; a < N; ++a)
        sum += (a % 2) ? 4 * f(bottom + a * h) : 2 * f(bottom + a * h);

    sum += f(bottom) + f(top);
    sum *= h / 3;
    return sum;
}

double V(double steps, double X, function<double(double, double, double, double)> const& f) {
    return simpson(-1, 1, steps, [X, steps, f](double x) {
        return simpson(-1, 1, steps, [X, steps, x, f](double y) {
            return simpson(-1, 1, steps, [X, x, y, f](double z) { return f(x, y, z, X); });
        });
    });
}

int main() {
    int steps = 100;

    ofstream file1;
    ofstream file2;
    file1.open("V1.txt");
    file2.open("V2.txt");

    file1 << "X\tV\tMonopol" << endl;
    file2 << "X\tV\tDipol" << endl;

    for (int n = 0; n <= 80; n++) {
        double X = n * 0.1;
        file1 << X << "\t"
              <<  V(steps, X, [](double x, double y, double z, double X) {
                      if (X - x < 1e-10 && y < 1e-10 && z < 1e-10) return 0.;
                      return 1./sqrt((X - x) * (X - x) + y * y + z * z);
                  })
              << "\t"
              << (1./X) * V(steps, X, [](double x, double y, double z, double X) { return 1; })
              << endl;

        file2 << X << "\t"
              <<  V(steps, X, [](double x, double y, double z, double X) {
                      if (X - x < 1e-10 && y < 1e-10 && z < 1e-10) return 0.;
                      return x / sqrt((X - x) * (X - x) + y * y + z * z);
                  })
              << "\t"
              << 1./(X*X) * V(steps, X, [](double x, double y, double z, double X) { return x*x; })
              << endl;
    }
    file1.close();
    file2.close();

    return 0;
}
