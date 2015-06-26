#include <random>
#include <fstream>
#include <iostream>
#include <functional>
#include <math.h>

long double mc_int_over_ellipse(unsigned long long N, double a, double b,
                                const std::function<double(double, double)> &f) {
    static unsigned long long seed(12345678912345678912ull);
    static std::mt19937_64 gen1(seed >> 1);
    static std::mt19937_64 gen2(seed << 1);
    std::uniform_real_distribution<long double> dist1(0.0, a);
    std::uniform_real_distribution<long double> dist2(0.0, b);

    long double sum = 0;

    for (unsigned int i = 0; i < N; i++) {
        long double x = dist1(gen1);
        long double y = dist2(gen2);

        if (x * x / (a * a) + y * y / (b * b) < 1) sum += f(x, y);
    }

    return sum * a * b * 4. / N;
}

double rnd_walk_pi(size_t N, unsigned int N_0, double a,
                   const std::function<double(double, double)> &boundary = [](double x, double) { return x; }) {
    unsigned long long seed(12345678912345678912ull);
    std::mt19937_64 gen1(seed << 1);
    std::mt19937_64 gen2(seed >> 1);
    std::uniform_real_distribution<double> dist(-a, std::nextafter(a, a + 1));

    unsigned long sum = 0;

    double x(0), y(0), delta;
    for (unsigned int i = 0; i < N; ++i) {
        delta = dist(gen1);
        if (fabs(x + delta) > 1)
            x = boundary(x, delta);
        else
            x += delta;
        delta = dist(gen2);
        if (fabs(y + delta) > 1)
            y = boundary(y, delta);
        else
            y += delta;

        if (i > N_0 && x * x + y * y < 1) sum++;
        // if (i < 1000) std::cout << x << "\t" << y << std::endl;
    }

    return sum * 4. / (N - N_0);
}

int main() {
    // a,b,c,d
    std::ofstream file("mc_int.txt");
    file << "# N\tpi(approx)\tdelta\tellipsen\tf_int" << std::endl;
    file.precision(5);

    for (int k = 1; k < 7; k++) {
        unsigned long long N = pow(10, k);
        long double pi = mc_int_over_ellipse(N, 1, 1, [](double, double) { return 1.; });     // pi
        long double ell_1 = mc_int_over_ellipse(N, 1, 1, [](double, double) { return 1.; });  // 1*pi
        long double ell_2 = mc_int_over_ellipse(N, 2, 1, [](double, double) { return 1.; });  // 2*pi
        long double ell_3 = mc_int_over_ellipse(N, 3, 1, [](double, double) { return 1.; });  // 3*pi
        long double ell_4 = mc_int_over_ellipse(N, 2, 2, [](double, double) { return 1.; });  // 4*pi
        long double ell_5 = mc_int_over_ellipse(N, 5, 1, [](double, double) { return 1.; });  // 5*pi
        long double ell_6 = mc_int_over_ellipse(N, 2, 3, [](double, double) { return 1.; });  // 6*pi
        long double ell_7 = mc_int_over_ellipse(N, 7, 1, [](double, double) { return 1.; });  // 7*pi
        long double ell_8 = mc_int_over_ellipse(N, 2, 4, [](double, double) { return 1.; });  // 8*pi
        long double ell_9 = mc_int_over_ellipse(N, 3, 3, [](double, double) { return 1.; });  // 9*pi
        long double f_int =
            mc_int_over_ellipse(N, sqrt(2), 1, [](double x, double) { return exp(-x * x); });  // 2.993...

        file << N << "\t" << pi << "\t" << M_PI - pi << "\t" << ell_1 << "\t" << ell_2 << "\t" << ell_3 << "\t" << ell_4
             << "\t" << ell_5 << "\t" << ell_6 << "\t" << ell_7 << "\t" << ell_8 << "\t" << ell_9 << "\t" << f_int
             << std::endl;

        std::cout << "d) Value of the integral is I = " << f_int << std::endl;
    }

    file.close();

    // b histogram
    file.open("mc_int_hist.txt");
    file.precision(10);

    for (int i = 0; i < 1000; ++i)
        file << mc_int_over_ellipse(1000, 1, 1, [](double, double) { return 1.; }) << std::endl;

    // e
    file.close();
    file.open("mc_int_rndwalk.txt");
    file.precision(5);

    for (int i = 3; i < 7; ++i) {
        unsigned int N = pow(10, i);
        file << N << "\t" << M_PI - rnd_walk_pi(N, pow(10, 2), 0.01) << "\t"
             << M_PI - rnd_walk_pi(N, pow(10, 2), 0.01, [](double x, double delta) {
                           return (x > 0) ? 1 - (fabs(x + delta) - 1) : -1 + (fabs(x + delta) - 1);
                       }) << std::endl;
    }

    return 0;
}