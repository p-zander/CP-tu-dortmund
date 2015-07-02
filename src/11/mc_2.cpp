#include <math.h>
#include <omp.h>
#include <random>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

using mattype = std::vector<std::vector<bool>>;

std::mt19937_64 gen;
constexpr size_t N = 100;
const double J = 1;

mattype init_spins_random() {
    mattype spins(N, std::vector<bool>(N));

    std::bernoulli_distribution dist;

    for (auto it = spins.begin(); it != spins.end(); ++it) {
        for (size_t j = 0; j < N; ++j) {
            (*it)[j] = dist(gen);
        }
    }

    return spins;
}

mattype init_spins_const() {
    mattype spins(N, std::vector<bool>(N));

    for (auto it = spins.begin(); it != spins.end(); ++it) {
        for (size_t j = 0; j < N; ++j) {
            (*it)[j] = 0;
        }
    }

    return spins;
}

int H(const mattype &spins, size_t x, size_t y) {  // over all 8 neighbors of element x,y (3x3)
    int sum = -8;

    for (int j = -1; j < 2; ++j) {
        for (int i = -1; i < 2; ++i) {
            if (i == j && i == 0) continue;
            sum += 2 * spins[(x + i) % N][(y + j) % N];
        }
    }
    return J * sum * ((spins[x][y]) ? 1 : -1);
}

double calc_E(const mattype &spins) {
    double E = 0;
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            E += H(spins, i, j);
        }
    }
    return E / (N * N);
}

double calc_M(const mattype &spins) {
    double M = -static_cast<double>(N * N);
    for (auto it = spins.begin(); it != spins.end(); ++it) {
        for (size_t j = 0; j < N; ++j) {
            M += 2 * (*it)[j];
        }
    }
    return M / (N * N);
}

double calc_Ck_BT2(const mattype &spins) {
    double sum(0), sum2(0);

    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            sum += H(spins, i, j);
            sum2 += pow(H(spins, i, j), 2);
        }
    }

    return sum2 / pow(N, 4) - pow(sum / (N * N), 2);
}

void simulate(mattype &spins, unsigned int steps, double k_BT) {
    std::uniform_int_distribution<size_t> x_dist(0, N - 1);
    std::uniform_int_distribution<size_t> y_dist(0, N - 1);

    std::uniform_real_distribution<double> dist(0., std::nextafter(1., 2.));

    for (unsigned int k = 0; k < steps; ++k) {
        size_t x = x_dist(gen);
        size_t y = y_dist(gen);

        int delta_E = 2 * H(spins, x, y);

        if (delta_E < 0)
            spins[x][y] = !spins[x][y];
        else if (dist(gen) < exp(-delta_E / k_BT))
            spins[x][y] = !spins[x][y];
    }
}

void print(const mattype &spins, std::string output_file) {
    std::ofstream file(output_file);

    for (auto it = spins.begin(); it != spins.end(); ++it) {
        for (size_t j = 0; j < (*it).size(); ++j) {
            file << (*it)[j] << " ";
        }
        file << std::endl;
    }
}

void a() {
    unsigned int steps = 1e6;
    double k_BT = 1.;

    mattype spins = init_spins_random();
    print(spins, "mc_2_init.txt");

    simulate(spins, steps, k_BT);
    print(spins, "mc_2_snap_1.txt");

    k_BT = 3.;

    spins = init_spins_random();
    simulate(spins, steps, k_BT);
    print(spins, "mc_2_snap_3.txt");
}

void b(unsigned int steps_max, int delta_steps, double k_BT, std::string output_file) {
    mattype spins_1 = init_spins_random();
    mattype spins_2 = init_spins_const();

    std::ofstream file(output_file);
    for (unsigned int i = 0; i < steps_max / delta_steps; ++i) {
        simulate(spins_1, delta_steps, k_BT);
        simulate(spins_2, delta_steps, k_BT);
        file << (i + 1) * delta_steps << "\t" << calc_E(spins_1) << "\t" << calc_E(spins_2) << std::endl;
    }
}

void c(unsigned int steps, int delta_steps, unsigned int steps_equi, double k_BT, std::string output_file) {
    mattype spins = init_spins_random();
    simulate(spins, steps_equi, k_BT);

    std::ofstream file(output_file);
    for (unsigned int i = 0; i < steps / delta_steps; ++i) {
        simulate(spins, delta_steps, k_BT);
        file << (i + 1) * delta_steps << "\t" << calc_E(spins) << "\t" << calc_M(spins) << std::endl;
    }
}

void d(double k_BT, std::ofstream &file) {
    mattype spins = init_spins_random();
    simulate(spins, 105e4, k_BT);

    std::stringstream buf;
    buf << k_BT << "\t" << abs(calc_M(spins)) << "\t" << k_BT * k_BT * calc_Ck_BT2(spins) << std::endl;

#pragma omp critical
    file << buf.str();
}

int main() {
    std::ofstream file_d("mc_2_d.txt");
#pragma omp parallel
    {
        if (omp_get_thread_num() == 0 % omp_get_num_threads()) a();
        if (omp_get_thread_num() == 1 % omp_get_num_threads()) {
            b(50000, 500, 1., "mc_2_b_1_00.txt");
            b(50000, 500, 2.25, "mc_2_b_2_25.txt");
            b(50000, 500, 3., "mc_2_b_3_00.txt");
        }
        if (omp_get_thread_num() == 2 % omp_get_num_threads()) c(1e6, 1e4, 50000, 1., "mc_2_c_1_00.txt");
        if (omp_get_thread_num() == 3 % omp_get_num_threads()) c(1e6, 1e4, 50000, 2.25, "mc_2_c_2_25.txt");
        if (omp_get_thread_num() == 4 % omp_get_num_threads()) c(1e6, 1e4, 50000, 3., "mc_2_c_3_00.txt");
        if (omp_get_thread_num() == 5 % omp_get_num_threads()) d(1.00, file_d);
        if (omp_get_thread_num() == 6 % omp_get_num_threads()) d(2.00, file_d);
        if (omp_get_thread_num() == 7 % omp_get_num_threads()) d(2.25, file_d);
        if (omp_get_thread_num() == 8 % omp_get_num_threads()) d(3.00, file_d);
        if (omp_get_thread_num() == 9 % omp_get_num_threads()) d(3.50, file_d);
        if (omp_get_thread_num() == 10 % omp_get_num_threads()) d(4.00, file_d);
        if (omp_get_thread_num() == 11 % omp_get_num_threads()) d(4.50, file_d);
    }

    return 0;
}
