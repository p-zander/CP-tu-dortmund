#include <math.h>
#include <random>
#include <fstream>

int main() {
    unsigned int N = 1e4;

    std::mt19937 gen;
    std::uniform_real_distribution<double> dist(0., std::nextafter(1., 2.));

    double k_BT = 1.;

    std::ofstream file("mc_1.txt");
    file.precision(10);

    for (int j = 1; j < 1000; ++j) {
        double H = 0.01 * j;

        short s = 1;
        int sum = 0;

        for (unsigned int i = 0; i < N; ++i) {
            if (s < 0)
                s *= -1;
            else if (dist(gen) < exp(-2 * H / k_BT))
                s *= -1;

            sum += s;
        }

        file << H << "\t" << sum * 1. / N << "\n";
    }
    file.close();
}
