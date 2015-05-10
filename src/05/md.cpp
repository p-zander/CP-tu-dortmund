#include <Eigen/Core>
#include <boost/random.hpp>
#include <boost/python.hpp>

#include <Python.h>

#include <sstream>
#include <iostream>
#include <utility>
#include <cmath>
#include <functional>
#include <algorithm>
#include <vector>
#include <string>

#include "../Utils/utils.h"

using std::pair;
using std::function;
using std::vector;
using std::cout;
using std::endl;
using std::make_pair;
using std::tie;
using Eigen::Vector2d;

double K_B = 1.3806488e-23;

Vector2d F(const Vector2d &x, double t) { return 24 * (pow(x.squaredNorm(), 3) - 2) / pow(x.squaredNorm(), 7) * x; }

// double F_norm(const Vector2d& x, double t) {
//     return -24./x.norm() * (2./pow(x.squaredNorm(), 6) -
//     1./pow(x.squaredNorm(), 3));
// }

Vector2d F(const vector<Vector2d> &particles, const Vector2d &x, double t, double L, double cutoff = -1) {
    if (cutoff <= 0) cutoff = L / 2;

    Vector2d F_eff(0, 0);
    for (const Vector2d &y : particles) {
        if (x == y) continue;
        for (int i = -1; i < 2; i++) {
            for (int j = -1; j < 2; j++) {
                Vector2d y_new = y + Vector2d(i * L, j * L);
                if ((x - y_new).norm() > cutoff) continue;
                F_eff += F(x - y_new, t);
            }
        }
    }
    return F_eff;
}

pair<Vector2d, Vector2d> Verlet_step_one(Vector2d x_0, Vector2d v_0, double m, double t_0, double h,
                                         const function<Vector2d(Vector2d, double)> &F) {
    Vector2d x_1 = x_0 + h * v_0 + h * h / (2 * m) * F(x_0, t_0);
    Vector2d v_1 = v_0 + h / (2 * m) * (F(x_0, t_0) + F(x_1, t_0 + h));
    return std::make_pair(x_1, v_1);
}

pair<vector<Vector2d>, vector<Vector2d>>
Verlet_step(const vector<Vector2d> &loc_0, const vector<Vector2d> &vel_0, double m, double t_0, double h,
            const function<Vector2d(vector<Vector2d>, Vector2d x, double)> &F) {
    assert(loc_0.size() == vel_0.size());
    vector<Vector2d> loc_1, vel_1;

    for (int i = 0; i < loc_0.size(); ++i) {
        Vector2d x, v;
        tie(x, v) = Verlet_step_one(loc_0[i], vel_0[i], m, t_0, h,
                                    [&loc_0, &F](Vector2d x, double t) { return F(loc_0, x, t); });
        loc_1.push_back(x);
        vel_1.push_back(v);
    }
    return make_pair(loc_1, vel_1);
}
// vector<pair<Vector2d, Vector2d>>
// Verlet(vector<Vector2d> x_0, vector<Vector2d> v_0, double m, double t_0,
//         double t_N, double h, const function<Vector2d(Vector2d, double)>& F) {
//     for (int i = 0; i <= t_N; i++) {
//     }
// }

pair<vector<Vector2d>, vector<Vector2d>> initialize(double L, double T_0, double m) {
    namespace rand = boost::random;

    rand::mt19937 generator;                                // Mersenne Twister Generator
    rand::uniform_real_distribution<> normdist(-1.0, 1.0);  // Distribution
    rand::variate_generator<rand::mt19937 &, rand::uniform_real_distribution<>> norm_rnd(generator, normdist);
    // Combination of distribution and generator

    vector<Vector2d> loc, vel;
    Vector2d v_mean(0, 0);
    double v_squared_mean = 0;
    int i = 0;

    for (int n = 0; n < 4; ++n) {
        for (int m = 0; m < 4; ++m) {
            loc.push_back(L / 8 * Vector2d(1 + 2 * n, 1 + 2 * m));

            Vector2d v_tmp = Vector2d(norm_rnd(), norm_rnd());
            vel.push_back(v_tmp);

            v_mean += v_tmp;
            v_squared_mean += v_tmp.squaredNorm();
            i++;
        }
    }
    v_mean /= i;
    v_squared_mean /= i;

    double scale = 10;  // 3 * K_B * T_0 / m;

    for (Vector2d &v : vel) v = (v - v_mean) * sqrt(scale / v_squared_mean);
    return make_pair(loc, vel);
}

int main(int argc, char const *argv[]) {
    int L = 8;
    double cutoff = L / 2;
    double m = 1;

    double h = 0.01;
    // double t_max = 10;

    vector<Vector2d> locations_0, locations_1;
    vector<Vector2d> velocities_0, velocities_1;

    tie(locations_0, velocities_0) = initialize(L, 1, m);

    tie(locations_1, velocities_1) =
        Verlet_step(locations_0, velocities_0, m, 0, h, [L, cutoff](const vector<Vector2d> &loc, const Vector2d &x,
                                                                    double t) { return F(loc, x, t, L, cutoff); });

    // prepare for plotting
    vector<double> x_0, y_0, v_x_0, v_y_0, F_x_0, F_y_0;
    vector<double> x_1, y_1, v_x_1, v_y_1, F_x_1, F_y_1;

    for (const Vector2d &r : locations_0) {
        x_0.push_back(r[0]);
        y_0.push_back(r[1]);

        Vector2d force = F(locations_0, r, 1, L, cutoff);
        F_x_0.push_back(force[0]);  // to test if the force is zero in the beginning
        F_y_0.push_back(force[1]);
    }
    for (const Vector2d &v : velocities_0) {
        v_x_0.push_back(v[0]);
        v_y_0.push_back(v[1]);
    }
    for (const Vector2d &r : locations_1) {
        x_1.push_back(r[0]);
        y_1.push_back(r[1]);

        Vector2d force = F(locations_1, r, 1, L, cutoff);
        F_x_1.push_back(force[0]);
        F_y_1.push_back(force[1]);
    }
    for (const Vector2d &v : velocities_1) {
        v_x_1.push_back(v[0]);
        v_y_1.push_back(v[1]);
    }

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

        utils::to_numpy(make_pair(x_0, y_0), make_pair("x_0", "y_0"), &global, &global);
        utils::to_numpy(make_pair(v_x_0, v_y_0), make_pair("v_x_0", "v_y_0"), &global, &global);
        utils::to_numpy(make_pair(F_x_0, F_y_0), make_pair("F_x_0", "F_y_0"), &global, &global);
        utils::to_numpy(make_pair(x_1, y_1), make_pair("x_1", "y_1"), &global, &global);
        utils::to_numpy(make_pair(v_x_1, v_y_1), make_pair("v_x_1", "v_y_1"), &global, &global);
        utils::to_numpy(make_pair(F_x_1, F_y_1), make_pair("F_x_1", "F_y_1"), &global, &global);

        // Launch some function in Python.
        py::exec_file("plots.py", global, global);
    } catch (const py::error_already_set &) {
        cout << utils::extractPythonException() << endl;
        exit(EXIT_FAILURE);
    }

    return 0;
}
