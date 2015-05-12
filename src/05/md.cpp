#include <Eigen/Core>
#include <boost/random.hpp>
#include <boost/python.hpp>

#include <Python.h>

#include <iostream>
#include <utility>
#include <functional>
#include <vector>

#include "../Utils/utils.h"

using std::pair;
using std::function;
using std::vector;
using std::cout;
using std::endl;
using std::make_pair;
using std::tie;
using std::bind;
using Eigen::Vector2d;

const double M = 1;
const double L = 10;
const unsigned int N = static_cast<unsigned int>(pow(3, 2));

inline Vector2d F_a_one(const Vector2d &x, double t) {
    // return 24 * (pow(x.squaredNorm(), 3) - 2) / pow(x.squaredNorm(), 7) * x;
    return 48 * x / x.squaredNorm() * (1. / pow(x.squaredNorm(), 6) - 1. / 2 * 1. / pow(x.squaredNorm(), 3));
}

Vector2d F_a(const vector<Vector2d> &particles, const Vector2d &x, double t, double cutoff) {
    Vector2d F_eff(0, 0);
    for (const Vector2d &y : particles) {
        if (x == y) continue;
        for (int i = -1; i < 2; i++) {
            for (int j = -1; j < 2; j++) {
                Vector2d y_new = y + Vector2d(i * L, j * L);
                if ((x - y_new).norm() > cutoff) continue;
                F_eff += F_a_one(y_new - x, t);
            }
        }
    }
    return F_eff;
}

pair<Vector2d, Vector2d> Verlet_step_one(Vector2d x_0, Vector2d v_0, double t_0, double h,
                                         const function<Vector2d(Vector2d, double)> &F) {
    Vector2d x_1 = x_0 + h * v_0 + h * h / (2 * M) * F(x_0, t_0);
    Vector2d v_1 = v_0 + h / (2 * M) * (F(x_0, t_0) + F(x_1, t_0 + h));
    return std::make_pair(x_1, v_1);
}

pair<vector<Vector2d>, vector<Vector2d>> Verlet_step(const vector<Vector2d> &loc_0, const vector<Vector2d> &vel_0,
                                                     double t_0, double h,
                                                     const function<Vector2d(vector<Vector2d>, Vector2d, double)> &F) {
    if (loc_0.size() != vel_0.size()) {
        std::cerr << "Verlet_step: Vectors loc_0 and vel_0 should have equal size! EXIT" << endl;
        exit(EXIT_FAILURE);
    }
    vector<Vector2d> loc_1, vel_1;

    for (size_t i = 0; i < loc_0.size(); ++i) {
        Vector2d x, v;
        tie(x, v) = Verlet_step_one(loc_0[i], vel_0[i], t_0, h, bind(F, loc_0, _1, _2));
        loc_1.push_back(x);
        vel_1.push_back(v);
    }
    return make_pair(loc_1, vel_1);
}

// vector<pair<Vector2d, Vector2d>>
// Verlet(vector<Vector2d> x_0, vector<Vector2d> v_0, double t_0,
//         double t_N, double h, const function<Vector2d(Vector2d, double)>& F) {
//     for (int i = 0; i <= t_N; i++) {
//     }
// }

pair<vector<Vector2d>, vector<Vector2d>> initialize(double T_0) {
    int ppa = static_cast<int>(sqrt(N));
    double dist = L / (2 * ppa);

    if (pow(ppa, 2) != N) {
        std::cerr << "initialize: sqrt(N) should be an integer! EXIT" << endl;
        exit(EXIT_FAILURE);
    }

    namespace rand = boost::random;

    rand::mt19937 generator;                                // Mersenne Twister Generator
    rand::uniform_real_distribution<> normdist(-1.0, 1.0);  // Distribution
    rand::variate_generator<rand::mt19937 &, rand::uniform_real_distribution<>> norm_rnd(generator, normdist);
    // Combination of distribution and generator

    vector<Vector2d> loc, vel;
    Vector2d v_mean(0, 0);
    double v_squared_mean = 0;

    for (int k = 0; k < ppa; ++k) {
        for (int l = 0; l < ppa; ++l) {
            loc.push_back(Vector2d(dist * (1 + 2 * k), dist * (1 + 2 * l)));

            Vector2d v_tmp = Vector2d(norm_rnd(), norm_rnd());
            vel.push_back(v_tmp);

            v_mean += v_tmp;
            v_squared_mean += v_tmp.squaredNorm();
        }
    }
    v_mean /= N;
    v_squared_mean /= N;

    double scale = 2 * T_0 / M;
    for (Vector2d &v : vel) v = (v - v_mean) * sqrt(scale / v_squared_mean);

    v_squared_mean = 0;
    for (const Vector2d &v : vel) v_squared_mean += v.squaredNorm();
    v_squared_mean /= N;

    cout << "temperature: " << .5 * M *v_squared_mean << endl;

    return make_pair(loc, vel);
}

int main(int argc, char const *argv[]) {
    double cutoff = L / 2 + 1;
    double h = 0.01;

    double t_0 = 0;
    // double t_max = 10;

    vector<Vector2d> locations_0, locations_1;
    vector<Vector2d> velocities_0, velocities_1;

    tie(locations_0, velocities_0) = initialize(0.5 /*= T_0*/);
    tie(locations_1, velocities_1) = Verlet_step(locations_0, velocities_0, h, h, bind(F_a, _1, _2, _3, cutoff));

    // prepare for plotting
    vector<Vector2d> F_0, F_1;

    for (const Vector2d &r : locations_0)
        F_0.push_back(F_a(locations_0, r, t_0, cutoff)); // to test if the force is zero in the beginning
    for (const Vector2d &r : locations_1)
        F_1.push_back(F_a(locations_1, r, t_0, cutoff));

    namespace py = boost::python;

    try {
        Py_Initialize();

        // Retrieve the main module's namespace
        py::object global(py::import("__main__").attr("__dict__"));

        py::exec("print 'Hello from Python!' \n"
                 "import numpy as np \n"
                 "print 'importing matplotlib...' \n"
                 "from matplotlib import pyplot as plt \n"
                 "print 'done' \n",
                 global, global);

        global["L"] = L;
        global["cutoff"] = cutoff;
        global["velocity_scale"] = 12;
        global["force_scale"] = 1e-3;

        utils::to_numpy(locations_0, "r_0", &global, &global);
        utils::to_numpy(velocities_0, "v_0", &global, &global);
        utils::to_numpy(F_0, "F_0", &global, &global);
        utils::to_numpy(locations_1, "r_1", &global, &global);
        utils::to_numpy(velocities_1, "v_1", &global, &global);
        utils::to_numpy(F_1, "F_1", &global, &global);

        // Launch some function in Python.
        py::exec_file("plots.py", global, global);
    } catch (const py::error_already_set &) {
        cout << utils::extractPythonException() << endl;
        exit(EXIT_FAILURE);
    }

    return 0;
}
