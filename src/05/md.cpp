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

using std::cout;
using std::endl;
using std::array;
using std::function;
using Eigen::Vector2d;

constexpr double M = 1;
constexpr double L = 8;
constexpr double cutoff = L / 2;
constexpr unsigned int N = 16;

Vector2d F_LJ(const Vector2d &x, double t) {
    double x2i = 1. / x.squaredNorm();
    double x6i = pow(x2i, 3);
    return 48 * x * x2i * x6i * (x6i - 0.5);
}

double V_LJ(const Vector2d &x, double t) {
    double x6i = 1. / pow(x.squaredNorm(), 3);
    return 4 * x6i * (x6i - 1);
}

std::pair<array<Vector2d, N>, double> forces(const array<Vector2d, N> &loc, double t) {
    array<Vector2d, N> F{{Vector2d(0, 0)}};
    double V = 0;

    for (size_t i = 0; i < N - 1; ++i) {
        for (size_t j = i + 1; j < N; ++j) {
            Vector2d r = loc[i] - loc[j];
            for (int k = -1; k < 2; ++k) {
                for (int m = -1; m < 2; ++m) {
                    Vector2d r_new = r + Vector2d(k * L, m * L);
                    if (r_new.norm() > cutoff) continue;
                    Vector2d force = F_LJ(r_new, t);
                    F[j] += force;
                    F[i] -= force;
                    V += V_LJ(r, t);
                }
            }
        }
    }

    return make_pair(F, V);
}

std::tuple<array<Vector2d, N>, array<Vector2d, N>, double>
Verlet_step(const array<Vector2d, N> &loc_0, const array<Vector2d, N> &vel_0, double t_0, double h) {
    if (loc_0.size() != vel_0.size()) {
        std::cerr << "Verlet_step: Vectors loc_0 and vel_0 should have equal size! EXIT" << endl;
        exit(EXIT_FAILURE);
    }
    array<Vector2d, N> loc_1, vel_1;

    array<Vector2d, N> F_0;
    std::tie(F_0, std::ignore) = forces(loc_0, t_0);

    for (size_t i = 0; i < N; ++i) {
        loc_1[i] = (loc_0[i] + h * vel_0[i] + h * h * F_0[i] / (2 * M))
                       .unaryExpr([](double xi) { return xi - L * floor(xi / L); });
    }
    double V_1;
    array<Vector2d, N> F_1;

    std::tie(F_1, V_1) = forces(loc_1, t_0);
    for (size_t i = 0; i < N; ++i) {
        vel_1[i] = vel_0[i] + h * (F_0[i] + F_1[i]) / (2 * M);
    }
    return std::make_tuple(loc_1, vel_1, V_1);
}

double Temp(const array<Vector2d, N> &vel) {
    double v_squared_mean = 0;
    for (const Vector2d &v : vel) v_squared_mean += v.squaredNorm();
    return M / 2 * v_squared_mean / N;
}

Vector2d CM_vel(const array<Vector2d, N> &vel) {
    return 1. / N * std::accumulate(vel.begin(), vel.end(), Vector2d(0, 0));
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
    double temp = Temp(vel);
    // cout << "Temp without scaling: " << temp << endl;
    double scale = sqrt(T_0 / temp);
    for (Vector2d &v : vel) v *= scale;

    Vector2d cm(CM_vel(vel));
    cout << "Init temperature: " << Temp(vel) << "\nInit CM velocity: ";

    cout << (cm.norm() < 1e-16) ? "0.0" : cm.norm() << endl;

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
    const double h = 0.001;

    const double t_0 = 0;
    const double t_max = 1e3;

    const double T_0 = 0.01;

    const size_t steps = static_cast<size_t>((t_max - t_0) / h + 0.5);
    const size_t samples = 10000;
    if (samples > steps) {
        std::cerr << "Number of samples should be less or equal number of steps! " << endl;
        exit(EXIT_FAILURE);
    }
    const size_t sbs = static_cast<size_t>(steps * 1. / samples);  // steps between samples

    typedef float dtype;
    boost::multi_array<dtype, 3> loc_sample(boost::extents[samples][N][2]);
    boost::multi_array<dtype, 3> vel_sample(boost::extents[samples][N][2]);
    boost::multi_array<dtype, 1> E_kin_sample(boost::extents[samples]);
    boost::multi_array<dtype, 1> E_pot_sample(boost::extents[samples]);
    boost::multi_array<dtype, 1> cm_vel_sample(boost::extents[samples]);

    array<array<Vector2d, N>, 2> locations, velocities;
    std::tie(locations[0], velocities[0]) = initialize(T_0);
    add_to_sample(locations[0], loc_sample, 0);
    add_to_sample(velocities[0], vel_sample, 0);

    cout << "Verlet in progress:  " << std::flush;
    double V;
    long at = 1;
    for (size_t t = 1; t < steps; t++) {
        std::tie(locations[t % 2], velocities[t % 2], V) =
            Verlet_step(locations[(t + 1) % 2], velocities[(t + 1) % 2], t_0 + t * h, h);
        if (t % sbs == 0) {
            printf((10 * t < steps) ? "%2.1f%%\b\b\b\b" : "%2.1f%%\b\b\b\b\b", static_cast<float>(100 * t / steps));
            E_kin_sample[at] = static_cast<float>(N * Temp(velocities[t % 2]));
            E_pot_sample[at] = static_cast<float>(V);
            cm_vel_sample[at] = static_cast<float>(CM_vel(velocities[t % 2]).norm());
            add_to_sample(locations[t % 2], loc_sample, at);
            add_to_sample(velocities[t % 2], vel_sample, at);
            at++;
        }
    }
    printf("%3.1f%%\ndone\n", 100.);

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
                 "from matplotlib import animation\n"
                 "plt.switch_backend('AGG') \n"
                 "print 'plotting...' \n",
                 global, global);

        global["file_ext"] = ".png";

        // Import variables and vectors
        global["L"] = L;
        global["N"] = N;
        global["samples"] = samples;
        global["cutoff"] = cutoff;

        global["velocity_scale"] = 12;
        global["force_scale"] = 1e-3;

        py::exec("strides_1 = np.ndarray((samples, N, 2), np.float32).strides", global, global);
        py::tuple strides_1 = py::extract<py::tuple>(global["strides_1"]);
        py::exec("strides_2 = np.ndarray((samples), np.float32).strides", global, global);
        py::tuple strides_2 = py::extract<py::tuple>(global["strides_2"]);

        global["loc_sample"] = np::from_data(loc_sample.data(), np::dtype::get_builtin<dtype>(),
                                             py::make_tuple(samples, N, 2), strides_1, py::object());
        global["vel_sample"] = np::from_data(vel_sample.data(), np::dtype::get_builtin<dtype>(),
                                             py::make_tuple(samples, N, 2), strides_1, py::object());
        global["E_kin_sample"] = np::from_data(E_kin_sample.data(), np::dtype::get_builtin<dtype>(),
                                               py::make_tuple(samples), strides_2, py::object());
        global["E_pot_sample"] = np::from_data(E_pot_sample.data(), np::dtype::get_builtin<dtype>(),
                                               py::make_tuple(samples), strides_2, py::object());
        global["cm_vel_sample"] = np::from_data(cm_vel_sample.data(), np::dtype::get_builtin<dtype>(),
                                                py::make_tuple(samples), strides_2, py::object());

        // Launch some function in Python.
        py::exec_file("plots.py", global, global);
    } catch (const py::error_already_set &) {
        cout << extractPythonException() << endl;
        exit(EXIT_FAILURE);
    }

    return 0;
}
