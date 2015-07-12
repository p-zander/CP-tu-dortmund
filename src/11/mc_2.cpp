#include <math.h>
#include <omp.h>
#include <random>
#include <vector>
#include <string>
#include <iostream>
#include <functional>
#include <fstream>
#include <sstream>

using spintype = bool;  // in order to access individual bits
using systype = std::vector<std::vector<spintype>>;

enum e { EandM, justM };

const double J = 1;

inline short spin(const spintype b) { return (b ? +1 : -1); }
inline size_t boundary(const size_t x, const short i, const size_t size) {
    return ((x == 0 && i == -1) ? (size - 1) : ((x == (size - 1) && i == 1) ? 0 : (x + i)));
}

systype init_sys(const size_t size, const std::function<spintype()> &f) {
    systype sys(size, std::vector<spintype>(size));

    for (auto it_1 = sys.begin(); it_1 != sys.end(); ++it_1) {
        for (auto it_2 = (*it_1).begin(); it_2 != (*it_1).end(); ++it_2) *it_2 = f();
    }

    return sys;
}

inline systype init_sys_const(const size_t size, const short val = 1) {
    return init_sys(size, [val]() { return val; });
}

inline systype init_sys_rand(const size_t size, std::mt19937_64 &gen, const float p = 0.5) {
    std::bernoulli_distribution dist(p);
    return init_sys(size, [&dist, &gen]() { return dist(gen); });
}

inline int sum_4(const systype &sys, const size_t x, const size_t y) {
    // over all 4 nearest neighbors of element x,y (3x3)
    return spin(sys[boundary(x, 1, sys.size())][y]) + spin(sys[x][boundary(y, 1, sys.size())]) +
           spin(sys[boundary(x, -1, sys.size())][y]) + spin(sys[x][boundary(y, -1, sys.size())]);
}

inline int sum_2(const systype &sys, const size_t x, const size_t y) {
    return spin(sys[boundary(x, 1, sys.size())][y]) + spin(sys[x][boundary(y, 1, sys.size())]);
}

// actual monte carlo simulation steps
void simulate(systype &sys, std::mt19937_64 &gen, const unsigned int steps, const double k_BT) {
    const size_t size = sys.size();

    std::uniform_int_distribution<size_t> x_dist(0, size - 1);
    std::uniform_int_distribution<size_t> y_dist(0, size - 1);

    std::uniform_real_distribution<double> p_dist(0., std::nextafter(1., 2.));

    for (unsigned int k = 0; k < steps; ++k) {
        size_t x = x_dist(gen);
        size_t y = y_dist(gen);

        double delta_E = 2 * J * spin(sys[x][y]) * sum_4(sys, x, y);

        if (delta_E < 0 || p_dist(gen) < exp(-delta_E / k_BT)) sys[x][y] = ((sys[x][y]) ? 0 : 1);
    }
}

// calculation of observables
double calc_E(const systype &sys) {
    size_t size = sys.size();
    double E = 0;

    for (size_t i = 0; i < size; ++i) {
        for (size_t j = 0; j < size; ++j) E += -J * spin(sys[i][j]) * sum_2(sys, i, j);
    }

    return E / (size * size);
}

double calc_M(const systype &sys) {
    size_t size = sys.size();
    double M = 0;

    for (auto it_1 = sys.begin(); it_1 != sys.end(); ++it_1) {
        for (auto it_2 = (*it_1).begin(); it_2 != (*it_1).end(); ++it_2) M += spin(*it_2);
    }
    return M / (size * size);
}

// several assignments
template <typename T> void print(const T &sys, std::string output_file) {
    std::ofstream file(output_file);

    for (auto it_1 = sys.begin(); it_1 != sys.end(); ++it_1) {
        for (auto it_2 = (*it_1).begin(); it_2 != (*it_1).end(); ++it_2) file << (*it_2) << " ";
        file << std::endl;
    }
}

void a(const size_t size, const unsigned long steps, std::mt19937_64 gen) {
    double k_BT = 1.;
    systype sys = init_sys_rand(size, gen);

    print(sys, "mc_2_init.txt");

    simulate(sys, gen, steps, k_BT);
    print(sys, "mc_2_snap_1.txt");

    k_BT = 3.;

    sys = init_sys_rand(size, gen);
    simulate(sys, gen, steps, k_BT);
    print(sys, "mc_2_snap_3.txt");
}

void b(const size_t size, const unsigned long steps_max, const unsigned long delta_steps, const double k_BT,
       std::mt19937_64 gen, const std::string output_file) {
    systype sys_rand = init_sys_rand(size, gen);
    systype sys_const = init_sys_const(size);

    std::ofstream file(output_file);
    file << 0 << "\t" << calc_E(sys_rand) << "\t" << calc_E(sys_const) << std::endl;

    for (unsigned int i = 0; i < steps_max / delta_steps; ++i) {
        simulate(sys_rand, gen, delta_steps, k_BT);
        simulate(sys_const, gen, delta_steps, k_BT);
        file << (i + 1) * delta_steps << "\t" << calc_E(sys_rand) << "\t" << calc_E(sys_const) << std::endl;
    }
}

void cde(size_t size, const unsigned long steps, const unsigned int delta_steps, const unsigned int steps_equi,
         const double k_BT, std::mt19937_64 gen, const std::string output_file, bool just_M = EandM) {
    systype sys = init_sys_rand(size, gen);
    simulate(sys, gen, steps_equi, k_BT);

    std::ofstream file(output_file);

    if (just_M) { file << 0 << "\t" << calc_M(sys) << std::endl; } else {
        file << 0 << "\t" << calc_E(sys) << "\t" << calc_M(sys) << std::endl;
    }

    for (unsigned int i = 0; i < steps / delta_steps; ++i) {
        simulate(sys, gen, delta_steps, k_BT);
        if (just_M) { file << (i + 1) * delta_steps << "\t" << calc_M(sys) << std::endl; } else {
            file << (i + 1) * delta_steps << "\t" << calc_E(sys) << "\t" << calc_M(sys) << std::endl;
        }
    }
}

int main() {
    unsigned long seed = 1337;

    const size_t N1 = 100;
    const size_t N2 = 50;
    const size_t N3 = 25;
    const size_t N4 = 20;
    const size_t N5 = 10;
    const size_t N6 = 5;
    const unsigned long sweep1 = N1 * N1;
    const unsigned long sweep2 = N2 * N2;
    const unsigned long sweep3 = N3 * N3;
    const unsigned long sweep4 = N4 * N4;
    const unsigned long sweep5 = N5 * N5;
    const unsigned long sweep6 = N6 * N6;

#pragma omp parallel
    {   
        std::mt19937_64 gen(seed + omp_get_thread_num());
            
        a(N1, 100 * sweep1, gen);
        if (omp_get_thread_num() == 0 % omp_get_num_threads()) {
        b(N1, 2e3 * sweep1, 10 * sweep1, 1.00, gen, "mc_2_b_100.txt");
        b(N1, 2e3 * sweep1, 10 * sweep1, 2.25, gen, "mc_2_b_225.txt");
        b(N1, 2e3 * sweep1, 10 * sweep1, 3.00, gen, "mc_2_b_300.txt");
        }

        if (omp_get_thread_num() == 0 % omp_get_num_threads())
            cde(N1, 1e5 * sweep1, 100 * sweep1, 2e3 * sweep1, 1.00, gen, "mc_2_cde_100.txt", EandM);
        if (omp_get_thread_num() == 1 % omp_get_num_threads())
            cde(N1, 1e5 * sweep1, 100 * sweep1, 2e3 * sweep1, 1.50, gen, "mc_2_cde_150.txt", EandM);
        if (omp_get_thread_num() == 2 % omp_get_num_threads())
            cde(N1, 1e5 * sweep1, 100 * sweep1, 2e3 * sweep1, 2.00, gen, "mc_2_cde_200.txt", EandM);
        if (omp_get_thread_num() == 3 % omp_get_num_threads())
            cde(N1, 1e5 * sweep1, 100 * sweep1, 2e3 * sweep1, 2.25, gen, "mc_2_cde_225.txt", EandM);
        if (omp_get_thread_num() == 4 % omp_get_num_threads())
            cde(N1, 1e5 * sweep1, 100 * sweep1, 2e3 * sweep1, 2.29, gen, "mc_2_cde_229.txt", EandM);
        if (omp_get_thread_num() == 5 % omp_get_num_threads())
            cde(N1, 1e5 * sweep1, 100 * sweep1, 2e3 * sweep1, 2.50, gen, "mc_2_cde_250.txt", EandM);
        if (omp_get_thread_num() == 6 % omp_get_num_threads())
            cde(N1, 1e5 * sweep1, 100 * sweep1, 2e3 * sweep1, 3.00, gen, "mc_2_cde_300.txt", EandM);
        if (omp_get_thread_num() == 7 % omp_get_num_threads())
            cde(N1, 1e5 * sweep1, 100 * sweep1, 2e3 * sweep1, 3.50, gen, "mc_2_cde_350.txt", EandM);

        if (omp_get_thread_num() == 1 % omp_get_num_threads()) {
            cde(N2, 1e5 * sweep2, 100 * sweep2, 5e2 * sweep2, 1.00, gen, "mc_2_U1_100.txt", justM);
            cde(N2, 1e5 * sweep2, 100 * sweep2, 5e2 * sweep2, 1.50, gen, "mc_2_U1_150.txt", justM);
            cde(N2, 1e5 * sweep2, 100 * sweep2, 5e2 * sweep2, 2.00, gen, "mc_2_U1_200.txt", justM);
            cde(N2, 1e5 * sweep2, 100 * sweep2, 5e2 * sweep2, 2.25, gen, "mc_2_U1_225.txt", justM);
        }
        if (omp_get_thread_num() == 2 % omp_get_num_threads()) {
            cde(N2, 1e5 * sweep2, 100 * sweep2, 5e2 * sweep2, 2.29, gen, "mc_2_U1_229.txt", justM);
            cde(N2, 1e5 * sweep2, 100 * sweep2, 5e2 * sweep2, 2.50, gen, "mc_2_U1_250.txt", justM);
            cde(N2, 1e5 * sweep2, 100 * sweep2, 5e2 * sweep2, 3.00, gen, "mc_2_U1_300.txt", justM);
            cde(N2, 1e5 * sweep2, 100 * sweep2, 5e2 * sweep2, 3.50, gen, "mc_2_U1_350.txt", justM);
        }
        if (omp_get_thread_num() == 3 % omp_get_num_threads()) {
            cde(N3, 1e5 * sweep3, 100 * sweep3, 5e2 * sweep3, 1.00, gen, "mc_2_U2_100.txt", justM);
            cde(N3, 1e5 * sweep3, 100 * sweep3, 5e2 * sweep3, 1.50, gen, "mc_2_U2_150.txt", justM);
            cde(N3, 1e5 * sweep3, 100 * sweep3, 5e2 * sweep3, 2.00, gen, "mc_2_U2_200.txt", justM);
            cde(N3, 1e5 * sweep3, 100 * sweep3, 5e2 * sweep3, 2.25, gen, "mc_2_U2_225.txt", justM);
            cde(N3, 1e5 * sweep3, 100 * sweep3, 5e2 * sweep3, 2.29, gen, "mc_2_U2_229.txt", justM);
            cde(N3, 1e5 * sweep3, 100 * sweep3, 5e2 * sweep3, 2.50, gen, "mc_2_U2_250.txt", justM);
            cde(N3, 1e5 * sweep3, 100 * sweep3, 5e2 * sweep3, 3.00, gen, "mc_2_U2_300.txt", justM);
            cde(N3, 1e5 * sweep3, 100 * sweep3, 5e2 * sweep3, 3.50, gen, "mc_2_U2_350.txt", justM);
        }
        if (omp_get_thread_num() == 4 % omp_get_num_threads()) {
            cde(N4, 1e5 * sweep4, 100 * sweep4, 5e2 * sweep4, 1.00, gen, "mc_2_U3_100.txt", justM);
            cde(N4, 1e5 * sweep4, 100 * sweep4, 5e2 * sweep4, 1.50, gen, "mc_2_U3_150.txt", justM);
            cde(N4, 1e5 * sweep4, 100 * sweep4, 5e2 * sweep4, 2.00, gen, "mc_2_U3_200.txt", justM);
            cde(N4, 1e5 * sweep4, 100 * sweep4, 5e2 * sweep4, 2.25, gen, "mc_2_U3_225.txt", justM);
            cde(N4, 1e5 * sweep4, 100 * sweep4, 5e2 * sweep4, 2.29, gen, "mc_2_U3_229.txt", justM);
            cde(N4, 1e5 * sweep4, 100 * sweep4, 5e2 * sweep4, 2.50, gen, "mc_2_U3_250.txt", justM);
            cde(N4, 1e5 * sweep4, 100 * sweep4, 5e2 * sweep4, 3.00, gen, "mc_2_U3_300.txt", justM);
            cde(N4, 1e5 * sweep4, 100 * sweep4, 5e2 * sweep4, 3.50, gen, "mc_2_U3_350.txt", justM);
            cde(N5, 1e5 * sweep5, 100 * sweep5, 5e2 * sweep5, 1.00, gen, "mc_2_U4_100.txt", justM);
            cde(N5, 1e5 * sweep5, 100 * sweep5, 5e2 * sweep5, 1.50, gen, "mc_2_U4_150.txt", justM);
            cde(N5, 1e5 * sweep5, 100 * sweep5, 5e2 * sweep5, 2.00, gen, "mc_2_U4_200.txt", justM);
            cde(N5, 1e5 * sweep5, 100 * sweep5, 5e2 * sweep5, 2.25, gen, "mc_2_U4_225.txt", justM);
            cde(N5, 1e5 * sweep5, 100 * sweep5, 5e2 * sweep5, 2.29, gen, "mc_2_U4_229.txt", justM);
            cde(N5, 1e5 * sweep5, 100 * sweep5, 5e2 * sweep5, 2.50, gen, "mc_2_U4_250.txt", justM);
            cde(N5, 1e5 * sweep5, 100 * sweep5, 5e2 * sweep5, 3.00, gen, "mc_2_U4_300.txt", justM);
            cde(N5, 1e5 * sweep5, 100 * sweep5, 5e2 * sweep5, 3.50, gen, "mc_2_U4_350.txt", justM);
            cde(N6, 1e5 * sweep6, 100 * sweep6, 5e2 * sweep6, 1.00, gen, "mc_2_U5_100.txt", justM);
            cde(N6, 1e5 * sweep6, 100 * sweep6, 5e2 * sweep6, 1.50, gen, "mc_2_U5_150.txt", justM);
            cde(N6, 1e5 * sweep6, 100 * sweep6, 5e2 * sweep6, 2.00, gen, "mc_2_U5_200.txt", justM);
            cde(N6, 1e5 * sweep6, 100 * sweep6, 5e2 * sweep6, 2.25, gen, "mc_2_U5_225.txt", justM);
            cde(N6, 1e5 * sweep6, 100 * sweep6, 5e2 * sweep6, 2.29, gen, "mc_2_U5_229.txt", justM);
            cde(N6, 1e5 * sweep6, 100 * sweep6, 5e2 * sweep6, 2.50, gen, "mc_2_U5_250.txt", justM);
            cde(N6, 1e5 * sweep6, 100 * sweep6, 5e2 * sweep6, 3.00, gen, "mc_2_U5_300.txt", justM);
            cde(N6, 1e5 * sweep6, 100 * sweep6, 5e2 * sweep6, 3.50, gen, "mc_2_U5_350.txt", justM);
        }
    }
}
