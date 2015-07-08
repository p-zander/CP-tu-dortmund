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

std::mt19937_64 gen;
const size_t N = 100;
const unsigned long sweep = N * N;
const double J = 1;

inline short spin(spintype b) { return (b ? +1 : -1); }
inline size_t boundary(const size_t x, const short i) {
    return ((x == 0 && i == -1) ? (N - 1) : ((x == (N - 1) && i == 1) ? 0 : (x + i)));
}

systype init_sys(const std::function<double()> &gen, size_t size = N) {
    systype sys(size, std::vector<spintype>(size));

    for (auto it = sys.begin(); it != sys.end(); ++it) {
        for (size_t j = 0; j < size; ++j) {
            (*it)[j] = gen();
        }
    }

    return sys;
}

inline systype init_sys_const(short val = 1) {
    return init_sys([val]() { return val; });
}

inline systype init_sys_rand(float p = 0.5, size_t size = N) {
    std::bernoulli_distribution dist(p);
    return init_sys([&dist]() { return dist(gen); }, size);
}

inline int sum_4(const systype &sys, const size_t x,
                 const size_t y) {  // over all 4 nearest neighbors of element x,y (3x3)
    return spin(sys[boundary(x, 1)][y]) + spin(sys[x][boundary(y, 1)]) + spin(sys[boundary(x, -1)][y]) +
           spin(sys[x][boundary(y, -1)]);
}

inline int sum_2(const systype &sys, const size_t x, const size_t y) {
    return spin(sys[boundary(x, 1)][y]) + spin(sys[x][boundary(y, 1)]);
}


// actual monte carlo simulation steps
void simulate(systype &sys, /*Htype &H,*/ const unsigned int steps, const double k_BT) {
    static std::uniform_int_distribution<size_t> x_dist(0, N - 1);
    static std::uniform_int_distribution<size_t> y_dist(0, N - 1);

    static std::uniform_real_distribution<double> dist(0., std::nextafter(1., 2.));
    for (unsigned int k = 0; k < steps; ++k) {
        size_t x = x_dist(gen);
        size_t y = y_dist(gen);

        double delta_E = 2 * J * spin(sys[x][y]) * sum_4(sys, x, y); 

        if (delta_E < 0 || dist(gen) < exp(-delta_E / k_BT)) {
            sys[x][y] = ((sys[x][y]) ? 0 : 1);
        }
    }
}

// calculation of observables
double calc_E(const systype &sys) {
    double E = 0;
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            E += -J * spin(sys[i][j]) * sum_2(sys, i, j);
        }
    }
    return E / (N * N);
}

double calc_M(const systype &sys) {
    double M = 0;
    for (auto it_1 = sys.begin(); it_1 != sys.end(); ++it_1) {
        for (auto it_2 = (*it_1).begin(); it_2 != (*it_1).end(); ++it_2) {
            M += spin(*it_2);
        }
    }
    return M / (N * N);
}


// several assignments
template <typename T> void print(const T &sys, std::string output_file) {
    std::ofstream file(output_file);

    for (auto it_1 = sys.begin(); it_1 != sys.end(); ++it_1) {
        for (auto it_2 = (*it_1).begin(); it_2 != (*it_1).end(); ++it_2) {
            file << (*it_2) << " ";
        }
        file << std::endl;
    }
}

void a() {
    unsigned int steps = 1000 * sweep;
    double k_BT = 1.;
    systype sys = init_sys_rand();

    print(sys, "mc_2_init.txt");

    simulate(sys, steps, k_BT);
    print(sys, "mc_2_snap_1.txt");

    k_BT = 3.;

    sys = init_sys_rand();
    simulate(sys, steps, k_BT);
    print(sys, "mc_2_snap_3.txt");
}

void b(unsigned long steps_max, unsigned long delta_steps, double k_BT, std::string output_file) {
    systype sys_rand = init_sys_rand();
    systype sys_const = init_sys_const();

    std::ofstream file(output_file);
    file << 0 << "\t" << calc_E(sys_rand) << "\t" << calc_E(sys_const) << std::endl;
    for (unsigned int i = 0; i < steps_max / delta_steps; ++i) {
        simulate(sys_rand, delta_steps, k_BT);
        simulate(sys_const, delta_steps, k_BT);
        file << (i + 1) * delta_steps << "\t" << calc_E(sys_rand) << "\t" << calc_E(sys_const) << std::endl;
    }
}

void cde(unsigned long steps, unsigned int delta_steps, unsigned int steps_equi, double k_BT, std::string output_file, bool just_M = false, size_t size = N){
    systype sys = init_sys_rand(size);
    simulate(sys,steps_equi, k_BT);


    std::ofstream file(output_file);
    
    if(just_M){
        file << 0 << "\t" << calc_M(sys) << std::endl;
        for(unsigned int i = 0; i < steps / delta_steps; ++i) {
            simulate(sys, delta_steps, k_BT);
            file << (i + 1) * delta_steps << "\t" << calc_M(sys) << std::endl;
        }
    }
    else{
        file << 0 << "\t" << calc_E(sys) << "\t" << calc_M(sys) << std::endl;
        for(unsigned int i = 0; i < steps / delta_steps; ++i) {
            simulate(sys, delta_steps, k_BT);
            file << (i + 1) * delta_steps << "\t" << calc_E(sys) << "\t" << calc_M(sys) << std::endl;
        }
    }
}


int main() {
    std::ofstream file_d("mc_2_d.txt");
    const size_t N2 = 50;
    const size_t N3 = 25;
    const unsigned long sweep2 = N2 * N2;
    const unsigned long sweep3 = N3 * N3;
#pragma omp parallel
    {
        if (omp_get_thread_num() == 0 % omp_get_num_threads()) {
            a();
            b(2e3 * sweep, 10*sweep, 1.00, "mc_2_b_1_00.txt");
            b(2e3 * sweep, 10*sweep, 2.25, "mc_2_b_2_25.txt");
            b(2e3 * sweep, 10*sweep, 3.00, "mc_2_b_3_00.txt");
        }
        if (omp_get_thread_num() == 1 % omp_get_num_threads()) cde(1e5*sweep, 100*sweep, 1.5e7, 1.00, "mc_2_cde_100.txt");
        if (omp_get_thread_num() == 2 % omp_get_num_threads()) cde(1e5*sweep, 100*sweep, 1.5e7, 1.50, "mc_2_cde_150.txt");
        if (omp_get_thread_num() == 3 % omp_get_num_threads()) cde(1e5*sweep, 100*sweep, 1.5e7, 2.00, "mc_2_cde_200.txt");
        if (omp_get_thread_num() == 4 % omp_get_num_threads()) cde(1e5*sweep, 100*sweep, 1.5e7, 2.25, "mc_2_cde_225.txt");
        if (omp_get_thread_num() == 5 % omp_get_num_threads()) cde(1e5*sweep, 100*sweep, 1.5e7, 2.50, "mc_2_cde_250.txt");
        if (omp_get_thread_num() == 6 % omp_get_num_threads()) cde(1e5*sweep, 100*sweep, 1.5e7, 3.00, "mc_2_cde_300.txt");
        if (omp_get_thread_num() == 7 % omp_get_num_threads()) cde(1e5*sweep, 100*sweep, 1.5e7, 3.50, "mc_2_cde_350.txt");

        if (omp_get_thread_num() == 8 % omp_get_num_threads()) cde(1e5*sweep2, 100*sweep2, 1.3e6, 1.00, "mc_2_U1_100.txt", true, N2);
        if (omp_get_thread_num() == 9 % omp_get_num_threads()) cde(1e5*sweep2, 100*sweep2, 1.3e6, 1.50, "mc_2_U1_150.txt", true, N2);
        if (omp_get_thread_num() == 10 % omp_get_num_threads()) cde(1e5*sweep2, 100*sweep2, 1.3e6, 2.00, "mc_2_U1_200.txt", true, N2);
        if (omp_get_thread_num() == 11 % omp_get_num_threads()) cde(1e5*sweep2, 100*sweep2, 1.3e6, 2.25, "mc_2_U1_225.txt", true, N2);
        if (omp_get_thread_num() == 12 % omp_get_num_threads()) cde(1e5*sweep2, 100*sweep2, 1.3e6, 2.50, "mc_2_U1_250.txt", true, N2);
        if (omp_get_thread_num() == 13 % omp_get_num_threads()) cde(1e5*sweep2, 100*sweep2, 1.3e6, 3.00, "mc_2_U1_300.txt", true, N2);
        if (omp_get_thread_num() == 14 % omp_get_num_threads()) cde(1e5*sweep2, 100*sweep2, 1.3e6, 3.50, "mc_2_U1_350.txt", true, N2);

        if (omp_get_thread_num() == 15 % omp_get_num_threads()) cde(1e5*sweep3, 100*sweep3, 1.2e6, 1.00, "mc_2_U2_100.txt", true, N3);
        if (omp_get_thread_num() == 16 % omp_get_num_threads()) cde(1e5*sweep3, 100*sweep3, 1.2e6, 1.50, "mc_2_U2_150.txt", true, N3);
        if (omp_get_thread_num() == 17 % omp_get_num_threads()) cde(1e5*sweep3, 100*sweep3, 1.2e6, 2.00, "mc_2_U2_200.txt", true, N3);
        if (omp_get_thread_num() == 18 % omp_get_num_threads()) cde(1e5*sweep3, 100*sweep3, 1.2e6, 2.25, "mc_2_U2_225.txt", true, N3);
        if (omp_get_thread_num() == 19 % omp_get_num_threads()) cde(1e5*sweep3, 100*sweep3, 1.2e6, 2.50, "mc_2_U2_250.txt", true, N3);
        if (omp_get_thread_num() == 20 % omp_get_num_threads()) cde(1e5*sweep3, 100*sweep3, 1.2e6, 3.00, "mc_2_U2_300.txt", true, N3);
        if (omp_get_thread_num() == 21 % omp_get_num_threads()) cde(1e5*sweep3, 100*sweep3, 1.2e6, 3.50, "mc_2_U2_350.txt", true, N3);

    }
}
