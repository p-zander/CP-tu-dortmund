#include <cmath>
#include <fstream>
#include <functional>
#include <tuple>

int main() {
    auto lcg_init = [](int64_t r0, int64_t a0, int64_t c0, int64_t m0) {
        return [=]() mutable {
            r0 = (a0 * r0 + c0) % m0;
            return r0 * 1. / m0;
        };
    };

    auto lcg1 = lcg_init(1234, 16807, 0, 2147483647);
    auto lcg2 = lcg_init(1234, 16807, 0, 2147483647);
    auto lcg3 = lcg_init(1234, 16807, 0, 2147483647);
    auto lcg4 = lcg_init(1234, 16807, 0, 2147483647);

    auto bm = [&lcg1]() {
        double x1 = lcg1();
        double x2 = lcg1();
        double y1 = sqrt(-2 * log(x1)) * cos(2 * M_PI * x2);
        double y2 = sqrt(-2 * log(x1)) * sin(2 * M_PI * x2);
        return std::make_tuple(y1, y2);
    };

    auto zgs = [&lcg2](unsigned int n) {
        double sum = 0;
        for (unsigned int i = 0; i < n; ++i) sum += (lcg2() - 0.5);
        return sum / sqrt(n * 1. / 12);
    };

    auto inv = [&lcg3]() { return pow(lcg3(), 1.0 / 3); };

    std::function<double()> rueck;
    rueck = [&lcg4, &rueck]() {
        auto x1 = lcg4();
        auto x2 = lcg4();
        if (x2 * M_PI < (sin(x1 * M_PI) / 2))
            return M_PI * x1;
        else
            return rueck();
    };

    double x2, x1;

    std::ofstream file("rand.txt");
    file << "#bm\tzgs\trueck\tinv" << std::endl;

    for (int i = 0; i < 100000; i++) {
        std::tie(x1, x2) = bm();
        file << x1 << "\t" << zgs(5000) << "\t" << rueck() << "\t" << inv() << std::endl;
        file << x2 << "\t" << zgs(5000) << "\t" << rueck() << "\t" << inv() << std::endl;
    }
    file.close();

    return 0;
}
