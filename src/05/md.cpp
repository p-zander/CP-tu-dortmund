#include <Eigen/Core>

#include <boost/random.hpp>
#include <boost/multi_array.hpp>
#include <boost/python.hpp>
#include <boost/numpy.hpp>

#include <Python.h>

#include <iostream>
#include <functional>
#include <utility>
#include <array>
#include <cmath>
#include <string>

// #include "../Utils/utils.h"

using std::cout;
using std::endl;
using std::array;
using std::function;
using Eigen::Vector2d;
namespace ph = std::placeholders;

constexpr double M = 1;
constexpr double L = 8;
constexpr unsigned int N = 16;

inline Vector2d F_LJ_one(const Vector2d &x, double t) { return -24 * x / (2 * pow(x.norm(), 14) - pow(x.norm(), -8)); }

Vector2d F_LJ(const array<Vector2d, N> &particles, const Vector2d &x, double t, double cutoff) {
    Vector2d F_eff(0, 0);
    for (const Vector2d &y : particles) {
        if (x.isApprox(y)) continue;
        for (int i = -1; i < 2; i++) {
            for (int j = -1; j < 2; j++) {
                Vector2d y_new = y + Vector2d(i * L, j * L);
                if ((x - y_new).norm() > cutoff) continue;
                F_eff += F_LJ_one(y_new - x, t);
            }
        }
    }
    return F_eff;
}

std::pair<Vector2d, Vector2d> Verlet_step_one(Vector2d x_0, Vector2d v_0, double t_0, double h,
                                              const function<Vector2d(Vector2d, double)> &F) {
    Vector2d x_1 = x_0 + h * v_0 + h * h / (2 * M) * F(x_0, t_0);
    x_1 = x_1.unaryExpr([](double xi) { return xi - L * floor(xi / L); });
    Vector2d v_1 = v_0 + h / (2 * M) * (F(x_0, t_0) + F(x_1, t_0 + h));
    return std::make_pair(x_1, v_1);
}

std::pair<array<Vector2d, N>, array<Vector2d, N>>
Verlet_step(const array<Vector2d, N> &loc_0, const array<Vector2d, N> &vel_0, double t_0, double h,
            const function<Vector2d(array<Vector2d, N>, Vector2d, double)> &F) {
    if (loc_0.size() != vel_0.size()) {
        std::cerr << "Verlet_step: Vectors loc_0 and vel_0 should have equal size! EXIT" << endl;
        exit(EXIT_FAILURE);
    }
    array<Vector2d, N> loc_1, vel_1;

    function<Vector2d(Vector2d, double)> F_tmp = std::bind(F, loc_0, ph::_1, ph::_2);
    for (size_t i = 0; i < loc_0.size(); ++i) {
        Vector2d x, v;
        std::tie(x, v) = Verlet_step_one(loc_0[i], vel_0[i], t_0, h, F_tmp);
        loc_1[i] = x;
        vel_1[i] = v;
    }
    return std::make_pair(loc_1, vel_1);
}

double Temp(const array<Vector2d, N> &vel) {
    double v_squared_mean = 0;
    for (const Vector2d &v : vel) v_squared_mean += v.squaredNorm();
    return M / 2 * v_squared_mean / N;
}

std::pair<array<Vector2d, N>, array<Vector2d, N>> initialize(double T_0) {
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

    array<Vector2d, N> loc, vel;

    Vector2d v_mean(0, 0);

    size_t i = 0;
    for (int k = 0; k < ppa; ++k) {
        for (int l = 0; l < ppa; ++l) {
            loc[i] = Vector2d(dist * (1 + 2 * k), dist * (1 + 2 * l));

            Vector2d v_tmp = Vector2d(norm_rnd(), norm_rnd());
            vel[i] = v_tmp;

            v_mean += v_tmp;

            i++;
        }
    }
    v_mean /= N;

    for (Vector2d &v : vel) v -= v_mean;
    double scale = sqrt(T_0 / Temp(vel));
    for (Vector2d &v : vel) v *= scale;

    cout << "Init temperature: " << Temp(vel) << endl;

    return std::make_pair(loc, vel);
}

std::string extractPythonException() {
    using namespace boost::python;

    PyObject *exc, *val, *tb;
    PyErr_Fetch(&exc, &val, &tb);
    PyErr_NormalizeException(&exc, &val, &tb);
    handle<> hexc(exc), hval(allow_null(val)), htb(allow_null(tb));

    if (!hval) {
        return extract<std::string>(str(hexc));
    } else {
        object traceback(import("traceback"));
        object format_exception(traceback.attr("format_exception"));
        object formatted_list(format_exception(hexc, hval, htb));
        object formatted(str("").join(formatted_list));
        return extract<std::string>(formatted);
    }
}

template <typename T1, typename T2>
void add_to_sample(const array<Eigen::Matrix<T1, 2, 1>, N> &vec, boost::multi_array<T2, 3> &sample, long at) {
    long i = 0;
    for (const Eigen::Matrix<T1, 2, 1> &r : vec) {
        sample[at][i][0] = static_cast<T2>(r[0]);
        sample[at][i][1] = static_cast<T2>(r[1]);
        i++;
    }
}

int main() {
    const double cutoff = L / 2;
    const double h = 0.01;

    const double t_0 = 0;
    const double t_max = 1e2;

    const double T_0 = 0.1;

    const size_t steps = static_cast<size_t>((t_max - t_0) / h + 0.5);
    const size_t samples = 10000;
    if (samples > steps) {
        cout << "ERROR: number samples should be less or equal number of steps! " << endl;
        exit(EXIT_FAILURE);
    }
    const size_t sbs = static_cast<size_t>(steps * 1. / samples);  // steps between samples

    typedef float dtype;
    typedef boost::multi_array<dtype, 3> array_type;
    typedef array_type::index index;
    array_type loc_sample(boost::extents[samples][N][2]);
    array_type vel_sample(boost::extents[samples][N][2]);

    array<array<Vector2d, N>, 2> locations, velocities;
    std::tie(locations[0], velocities[0]) = initialize(T_0);
    add_to_sample(locations[0], loc_sample, 0);
    add_to_sample(velocities[0], vel_sample, 0);

    cout << "Verlet in progress:  ";
    long at = 1;
    for (size_t t = 1; t < steps; t += 1) {
        std::tie(locations[t % 2], velocities[t % 2]) =
            Verlet_step(locations[(t + 1) % 2], velocities[(t + 1) % 2], t_0 + t * h, h,
                        std::bind(F_LJ, ph::_1, ph::_2, ph::_3, cutoff));
        if (t % sbs == 0) {
            printf((10*t < steps)? "%2.1f%%\b\b\b\b" : "%2.1f%%\b\b\b\b\b", static_cast<float>(100*t/steps));
            add_to_sample(locations[t % 2], loc_sample, at);
            add_to_sample(velocities[t % 2], vel_sample, at);
            at++;
        }
    }
    printf("%3.1f%%\ndone\n", 100.);

    // prepare for plotting
    //     array<Vector2d, N> F_0;

    //     size_t i = 0;
    //     auto F = std::bind(F_LJ, locations.back(), ph::_1, t_0, cutoff);
    //     for (const Vector2d &r : locations.back())
    //         F_0[i++] = F(r);  // to test if the force is zero in the beginning

    namespace py = boost::python;
    namespace np = boost::numpy;

    try {
        Py_Initialize();
        np::initialize();

        // Retrieve the main module's namespace
        py::object global(py::import("__main__").attr("__dict__"));

        // Import neccessary modules
        py::exec("print 'Hello from Python!' \n"
                 "import numpy as np \n"
                 "print 'importing matplotlib...' \n"
                 "from matplotlib import pyplot as plt \n"
                 "plt.switch_backend('AGG') \n"
                 "print 'plotting...' \n"
                 , global, global);

        global["file_ext"] = ".png";

        // Import variables and vectors
        global["L"] = L;
        global["N"] = N;
        global["samples"] = samples;
        global["cutoff"] = cutoff;

        global["velocity_scale"] = 12;
        global["force_scale"] = 1e-3;

        py::exec("strides = np.ndarray((samples, N, 2), np.float32).strides", global, global);
        py::tuple strides = py::extract<py::tuple>(global["strides"]);

        global["loc_sample"] = np::from_data(loc_sample.data(), np::dtype::get_builtin<dtype>(),
                                             py::make_tuple(samples, N, 2), strides, py::object());
        global["vel_sample"] = np::from_data(vel_sample.data(), np::dtype::get_builtin<dtype>(),
                                             py::make_tuple(samples, N, 2), strides, py::object());

        py::exec("plt.plot(loc_sample[0,:,0], loc_sample[0,:,1], 'bo') \n"
                 "plt.xlim(0, L) \nplt.ylim(0, L) \n"
                 "plt.savefig('init' + file_ext) \n"
                 "plt.figure() \n"
                 "for i in range(16): \n"
                 "    plt.plot(loc_sample[:,i,0], loc_sample[:,i,1], ',') \n"
                 "plt.savefig('time' + file_ext) \n",
                 global, global);

        // Launch some function in Python.
        // py::exec_file("plots.py", global, global);
        py::exec("print 'done' \n", global, global);
    } catch (const py::error_already_set &) {
        cout << extractPythonException() << endl;
        exit(EXIT_FAILURE);
    }

    return 0;
}
