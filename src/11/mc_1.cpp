#include <math.h>
#include <random>
#include <fstream>

int main() {
    unsigned int N = 1e4;

    std::mt19937_64 gen;
    std::uniform_real_distribution<double> dist(0., std::nextafter(1., 2.));

    double k_BT = 1.;

    std::ofstream file("mc_1.txt");
    file.precision(10);

    #pragma omp parallel for ordered
    for (int j = 0; j < 80; ++j) {
        double H = 0.05 * j;

        short s = 1;
        int sum = 0;

        for (unsigned int i = 0; i < N; ++i) {
            if (s < 0 || dist(gen) < exp(-2 * H / k_BT))
                s *= -1;

            sum += s;
        }

        #pragma omp ordered
        file << H << "\t" << sum * 1. / N << "\n";
    }
    file.close();
}
