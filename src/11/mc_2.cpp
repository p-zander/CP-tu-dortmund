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
using Htype = std::vector<std::vector<double>>;

std::mt19937_64 gen;
const size_t N = 20;
const unsigned long sweep = N * N;
const double J = 1;

inline short spin(spintype b) { return (b ? +1 : -1); }
inline size_t boundary(const size_t x, const short i) {
    return ((x == 0 && i == -1) ? (N - 1) : ((x == (N - 1) && i == 1) ? 0 : (x + i)));
}

systype init_sys(const std::function<double()> &gen) {
    systype sys(N, std::vector<spintype>(N));

    for (auto it = sys.begin(); it != sys.end(); ++it) {
        for (size_t j = 0; j < N; ++j) {
            (*it)[j] = gen();
        }
    }

    return sys;
}

inline systype init_sys_const(short val = 1) {
    return init_sys([val]() { return val; });
}

inline systype init_sys_rand(float p = 0.5) {
    std::bernoulli_distribution dist(p);
    return init_sys([&dist]() { return dist(gen); });
}

inline int sum_4(const systype &sys, const size_t x,
                 const size_t y) {  // over all 4 nearest neighbors of element x,y (3x3)
    return spin(sys[boundary(x, 1)][y]) + spin(sys[x][boundary(y, 1)]) + spin(sys[boundary(x, -1)][y]) +
           spin(sys[x][boundary(y, -1)]);
}

inline int sum_2(const systype &sys, const size_t x, const size_t y) {
    return spin(sys[boundary(x, 1)][y]) + spin(sys[x][boundary(y, 1)]);
}

// storage of the matrix of all local H(i) should be more efficient (depending on the temperature)
Htype calc_H_1(const systype &sys) {  // doesn't need the boundary. may be more efficient
    Htype H(N, std::vector<double>(N));
    int sum;

    // four corners
    sum = spin(sys[1][0]) + spin(sys[0][1]) + spin(sys[N - 1][0]) + spin(sys[0][N - 1]);
    H[0][0] = -J * spin(sys[0][0]) * sum;

    sum = spin(sys[0][0]) + spin(sys[N - 1][1]) + spin(sys[N - 2][0]) + spin(sys[N - 1][N - 1]);
    H[N - 1][0] = -J * spin(sys[N - 1][0]) * sum;

    sum = spin(sys[1][N - 1]) + spin(sys[0][0]) + spin(sys[0][N - 2]) + spin(sys[N - 1][N - 1]);
    H[0][N - 1] = -J * spin(sys[0][N - 1]) * sum;

    sum = spin(sys[0][N - 1]) + spin(sys[N - 1][0]) + spin(sys[N - 2][N - 1]) + spin(sys[N - 1][N - 2]);
    H[N - 1][N - 1] = -J * spin(sys[N - 1][N - 1]) * sum;

    // four edges
    for (size_t i = 1; i < N - 1; ++i) {
        sum = spin(sys[1][i]) + spin(sys[0][i + 1]) + spin(sys[N - 1][i]) + spin(sys[0][i - 1]);
        H[0][i] = -J * spin(sys[0][i]) * sum;

        sum = spin(sys[i + 1][0]) + spin(sys[i][1]) + spin(sys[i - 1][0]) + spin(sys[i][N - 1]);
        H[i][0] = -J * spin(sys[i][0]) * sum;

        sum = spin(sys[0][i]) + spin(sys[N - 1][i + 1]) + spin(sys[N - 2][i]) + spin(sys[N - 1][i - 1]);
        H[N - 1][i] = -J * spin(sys[N - 1][i]) * sum;

        sum = spin(sys[i + 1][N - 1]) + spin(sys[i][0]) + spin(sys[i - 1][N - 1]) + spin(sys[i][N - 2]);
        H[i][N - 1] = -J * spin(sys[i][N - 1]) * sum;
    }

    // inner
    for (size_t i = 1; i < N - 1; ++i) {
        for (size_t j = 1; j < N - 1; ++j) {
            sum = spin(sys[i + 1][j]) + spin(sys[i][j + 1]) + spin(sys[i - 1][j]) + spin(sys[i][j - 1]);
            H[i][j] = -J * spin(sys[i][j]) * sum;
        }
    }

    return H;
}

Htype calc_H_2(const systype &sys) {  // easier to read but does need the boundary conditions
    Htype H(N, std::vector<double>(N));

    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            H[i][j] = -J * spin(sys[i][j]) * sum_4(sys, i, j);
        }
    }

    return H;
}

// update the matrix H when a spin is flipped
void update_H(Htype &H, const systype &sys, const size_t x, const size_t y) {
    size_t k;

    H[x][y] *= -1;
    k = boundary(x, 1);
    H[k][y] = -J * spin(sys[k][y]) * sum_4(sys, k, y);
    k = boundary(y, 1);
    H[x][k] = -J * spin(sys[x][k]) * sum_4(sys, x, k);
    k = boundary(x, -1);
    H[k][y] = -J * spin(sys[k][y]) * sum_4(sys, k, y);
    k = boundary(y, -1);
    H[x][k] = -J * spin(sys[x][k]) * sum_4(sys, x, k);
}

// actual monte carlo simulation steps
void simulate(systype &sys, Htype &H, const unsigned int steps, const double k_BT) {
    static std::uniform_int_distribution<size_t> x_dist(0, N - 1);
    static std::uniform_int_distribution<size_t> y_dist(0, N - 1);

    static std::uniform_real_distribution<double> dist(0., std::nextafter(1., 2.));

    for (unsigned int k = 0; k < steps; ++k) {
        size_t x = x_dist(gen);
        size_t y = y_dist(gen);

        double delta_E = -2 * H[x][y];
        // double delta_E = 2 * J * spin(sys[x][y]) * sum_4(sys, x, y);  // is actually still faster

        if (delta_E < 0 || dist(gen) < exp(-delta_E / k_BT)) {
            sys[x][y] = ((sys[x][y]) ? 0 : 1);
            update_H(H, sys, x, y);  // to much array accessing
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

double calc_Ck_BT2(const Htype &H) {
    double sum(0), sum2(0);

    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            sum += H[i][j];
            sum2 += pow(H[i][j], 2);
        }
    }

    return sum2 / pow(N, 2) - pow(sum / (N * N), 2);
}

double calc_U(const systype &sys) {
    double sum2 = 0;
    double sum4 = 0;
    for (auto it_1 = sys.begin(); it_1 != sys.end(); ++it_1) {
        for (auto it_2 = (*it_1).begin(); it_2 != (*it_1).end(); ++it_2) {
            sum2 += pow(spin(*it_2), 2);
            sum4 += pow(spin(*it_2), 4);
        }
    }
    return 1 - (sum4 / (N * N)) / (3 * pow(sum2 / (N * N), 2));
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
    Htype H = calc_H_1(sys);
    print(sys, "mc_2_init.txt");

    simulate(sys, H, steps, k_BT);
    print(sys, "mc_2_snap_1.txt");

    k_BT = 3.;

    sys = init_sys_rand(0);
    H = calc_H_1(sys);
    simulate(sys, H, steps, k_BT);
    print(sys, "mc_2_snap_3.txt");
}

void b(unsigned long steps_max, unsigned long delta_steps, double k_BT, std::string output_file) {
    systype sys_rand = init_sys_rand();
    Htype H_rand = calc_H_1(sys_rand);
    systype sys_const = init_sys_const();
    Htype H_const = calc_H_1(sys_const);

    std::ofstream file(output_file);
    file << 0 << "\t" << calc_E(sys_rand) << "\t" << calc_E(sys_const) << std::endl;
    for (unsigned int i = 0; i < steps_max / delta_steps; ++i) {
        simulate(sys_rand, H_rand, delta_steps, k_BT);
        simulate(sys_const, H_const, delta_steps, k_BT);
        file << (i + 1) * delta_steps << "\t" << calc_E(sys_rand) << "\t" << calc_E(sys_const) << std::endl;
    }
}

void c(unsigned long steps, unsigned int delta_steps, unsigned int steps_equi, double k_BT, std::string output_file) {
    systype sys = init_sys_rand();
    Htype H = calc_H_1(sys);
    simulate(sys, H, steps_equi, k_BT);

    std::ofstream file(output_file);
    file << 0 << "\t" << calc_E(sys) << "\t" << calc_M(sys) << std::endl;
    for (unsigned int i = 0; i < steps / delta_steps; ++i) {
        simulate(sys, H, delta_steps, k_BT);
        file << (i + 1) * delta_steps << "\t" << calc_E(sys) << "\t" << calc_M(sys) << std::endl;
    }
}

void d(double k_BT, unsigned long steps, std::ofstream &file) {
    systype sys = init_sys_rand();
    Htype H = calc_H_1(sys);
    simulate(sys, H, steps, k_BT);

    std::stringstream buf;
    buf << k_BT << "\t" << fabs(calc_M(sys)) << "\t" << k_BT * k_BT * calc_Ck_BT2(H) << "\t" << calc_U(sys)
        << std::endl;

#pragma omp critical
    file << buf.str();
}

int main() {
    std::ofstream file_d("mc_2_d.txt");
#pragma omp parallel
    {
        if (omp_get_thread_num() == 0 % omp_get_num_threads()) {
            a();
            b(1e2 * sweep, sweep, 1.00, "mc_2_b_1_00.txt");
            b(1e2 * sweep, sweep, 2.25, "mc_2_b_2_25.txt");
            b(1e2 * sweep, sweep, 3.00, "mc_2_b_3_00.txt");
            c(1e3 * sweep, sweep, 1e6, 1.00, "mc_2_c_1_00.txt");
            c(1e3 * sweep, sweep, 1e6, 2.25, "mc_2_c_2_25.txt");
            c(1e3 * sweep, sweep, 1e6, 3.00, "mc_2_c_3_00.txt");
        }
        if (omp_get_thread_num() == 1 % omp_get_num_threads()) d(1.00 /*2.10*/, 3e5 + 1e5 * sweep, file_d);
        if (omp_get_thread_num() == 2 % omp_get_num_threads()) d(2.00 /*2.15*/, 3e5 + 1e5 * sweep, file_d);
        if (omp_get_thread_num() == 3 % omp_get_num_threads()) d(2.20 /*2.20*/, 3e5 + 1e5 * sweep, file_d);
        if (omp_get_thread_num() == 4 % omp_get_num_threads()) d(2.30 /*2.25*/, 3e5 + 1e5 * sweep, file_d);
        if (omp_get_thread_num() == 5 % omp_get_num_threads()) d(2.50 /*2.30*/, 3e5 + 1e5 * sweep, file_d);
        if (omp_get_thread_num() == 6 % omp_get_num_threads()) d(3.00 /*2.35*/, 3e5 + 1e5 * sweep, file_d);
        if (omp_get_thread_num() == 7 % omp_get_num_threads()) d(4.00 /*2.40*/, 3e5 + 1e5 * sweep, file_d);
    }

    return 0;
}
