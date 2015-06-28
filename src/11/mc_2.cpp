#include <math.h>
#include <omp.h>
#include <random>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

using mattype = std::vector<std::vector<bool>>;

const double J = 1;

mattype init_spins_random(size_t N) {
    mattype spins(N, std::vector<bool>(N));

    std::mt19937 gen;
    std::bernoulli_distribution dist;

    for (auto it = spins.begin(); it != spins.end(); ++it) {
        for (size_t j = 0; j < N; ++j) {
            (*it)[j] = dist(gen);
        }
    }

    return spins;
}

mattype init_spins_const(size_t N) {
    mattype spins(N, std::vector<bool>(N));

    std::mt19937 gen;
    std::bernoulli_distribution dist;

    for (auto it = spins.begin(); it != spins.end(); ++it) {
        for (size_t j = 0; j < N; ++j) {
            (*it)[j] = 1;
        }
    }

    return spins;
}

int H(const mattype &spins, size_t x, size_t y) {  // over all 8 neighbors of element x,y (3x3)
    int sum = -8;
    size_t N = spins.size();

    for (int i = -1; i < 2; ++i) {
        for (int j = -1; j < 2; ++j) {
            if (i == j && i == 0) continue;
            sum += 2 * spins[(x + i) % N][(y + j) % N];
        }
    }
    return -J * sum * ((spins[x][y]) ? 1 : -1);
}

double calc_E(const mattype &spins) {
    double E = 0;
    size_t N = spins.size();
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            E += H(spins, i, j);
        }
    }
    return E / (N * N);
}

double calc_M(const mattype &spins) {
    size_t N = spins.size();
    double M = -static_cast<double>(N * N);
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            M += 2 * spins[i][j];
        }
    }
    return M / static_cast<double>(N * N);
}

void simulate(mattype &spins, unsigned int steps, double k_BT) {
    static std::mt19937_64 loc_gen;
    static std::uniform_int_distribution<size_t> x_dist(0, spins.size() - 1);
    static std::uniform_int_distribution<size_t> y_dist(0, spins.size() - 1);

    static std::mt19937 gen;
    static std::uniform_real_distribution<double> dist(0., std::nextafter(1., 2.));

    for (unsigned int k = 0; k < steps; ++k) {
        size_t x = x_dist(loc_gen);
        size_t y = y_dist(loc_gen);

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
    size_t N = 100;
    double k_BT = 1.;

    mattype spins = init_spins_random(N);

    simulate(spins, steps, k_BT);
    print(spins, "mc_2_snap_1.txt");

    k_BT = 3.;

    spins = init_spins_random(N);
    simulate(spins, steps, k_BT);
    print(spins, "mc_2_snap_3.txt");
}

void b(unsigned int steps_max, int delta_steps, double k_BT, std::string output_file) {
    size_t N = 100;

    mattype spins_1 = init_spins_random(N);
    mattype spins_2 = init_spins_const(N);

    std::ofstream file(output_file);
    for (unsigned int i = 0; i < steps_max / delta_steps; ++i) {
        simulate(spins_1, delta_steps, k_BT);
        simulate(spins_2, delta_steps, k_BT);
        file << (i + 1) * delta_steps << "\t" << calc_E(spins_1) << "\t" << calc_E(spins_2) << std::endl;
    }
}

void c(unsigned int steps, int delta_steps, unsigned int steps_equi, double k_BT, std::string output_file) {
    size_t N = 100;

    mattype spins = init_spins_random(N);
    simulate(spins, steps_equi, k_BT);

    std::ofstream file(output_file);
    for (unsigned int i = 0; i < steps / delta_steps; ++i) {
        simulate(spins, delta_steps, k_BT);
        file << (i + 1) * delta_steps << "\t" << calc_E(spins) << "\t" << calc_M(spins) << std::endl;
    }
}

int main() {
    a();
#pragma omp parallel
    {
        if (omp_get_thread_num() == 0 % omp_get_num_threads()) b(30000, 300, 1., "mc_2_b_1_00.txt");
        if (omp_get_thread_num() == 1 % omp_get_num_threads()) b(30000, 300, 2.25, "mc_2_b_2_25.txt");
        if (omp_get_thread_num() == 2 % omp_get_num_threads()) b(30000, 300, 3., "mc_2_b_3_00.txt");
    }

#pragma omp parallel
    {
        if (omp_get_thread_num() == 0 % omp_get_num_threads()) c(1e6, 1e4, 50000, 1., "mc_2_c_1_00.txt");
        if (omp_get_thread_num() == 1 % omp_get_num_threads()) c(1e6, 1e4, 50000, 2.25, "mc_2_c_2_25.txt");
        if (omp_get_thread_num() == 2 % omp_get_num_threads()) c(1e6, 1e4, 50000, 3., "mc_2_c_3_00.txt");
    }

    return 0;
}
