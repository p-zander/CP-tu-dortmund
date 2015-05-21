#include <Eigen/Core>

#include <boost/random.hpp>
#include <boost/multi_array.hpp>
#include <boost/python.hpp>
#include <boost/numpy.hpp>
#include <boost/lexical_cast.hpp>

#include <Python.h>

#include <omp.h>

#include <iostream>
#include <tuple>
#include <array>
#include <cmath>
#include <string>

using Eigen::Vector2d;

constexpr double M = 1;
constexpr double L = 8;
constexpr double cutoff = L / 2;
constexpr unsigned int N = 16;

constexpr size_t bins = 500;

const double t_aqui = 10;

// Lennard-Jones force
inline Vector2d F_LJ(const Vector2d &x, double t) {
    double x2i = 1. / x.squaredNorm();
    double x6i = pow(x2i, 3);
    return 48 * x * x2i * x6i * (x6i - 0.5);
}

// Lennard-Jones potential
double V_LJ(const std::array<Vector2d, N> &x, double t) {
    double V = 0;
    Vector2d r;
    for (size_t i = 0; i < N - 1; ++i) {
        for (size_t j = i + 1; j < N; ++j) {
            r = (x[j] - x[i]).unaryExpr([](double xi) { return xi - L * floor(xi / L + 0.5); });
            if (r.norm() >= cutoff) continue;
            double x6i = 1. / pow(r.squaredNorm(), 3);
            V += 4 * x6i * (x6i - 1);
        }
    }
    return V;
}

// temperature, i.e. kinetical energy / N
double Temp(const std::array<Vector2d, N> &vel) {
    double v_squared_mean = 0;
    for (const Vector2d &v : vel) v_squared_mean += v.squaredNorm();
    return M / 2 * v_squared_mean / N;
}

// center of mass velocity
Vector2d CM_vel(const std::array<Vector2d, N> &vel) {
    return 1. / N * std::accumulate(vel.begin(), vel.end(), Vector2d(0, 0));
}

// calculation of all the N force vectors
std::array<Vector2d, N> calcForces(const std::array<Vector2d, N> &loc, double t) {
    std::array<Vector2d, N> F{{Vector2d(0, 0)}};
    Vector2d r;
    Vector2d force{0,0};

    for (size_t i = 0; i < N - 1; ++i) {
        for (size_t j = i + 1; j < N; ++j) {
            r = (loc[i] - loc[j]).unaryExpr([](double xi) { return xi - L * floor(xi / L + 0.5); });
            if (r.norm() >= cutoff) continue;
            force = F_LJ(r, t);
            F[i] += force;
            F[j] -= force;
        }
    }

    return F;
}

// update bins for g(r)
void update_radial(const std::array<Vector2d, N> &x, std::array<float, bins> &g, size_t bins) {
    Vector2d r;
    for (size_t i = 0; i < N - 1; ++i) {
        for (size_t j = i + 1; j < N; ++j) {
            r = (x[j] - x[i]).unaryExpr([](double xi) { return xi - L * floor(xi / L + 0.5); });
            if (r.norm() >= cutoff) continue;
            size_t index = size_t(floor(bins * r.norm() / cutoff));
            if (index > bins) {
                std::cerr << "ERROR: index for g is out of range. check computation! EXIT " << std::endl;
                exit(EXIT_FAILURE);
            }
            g[index] += 2;
        }
    }
}

// actual verlet step using given locations, velocities and forces from the next step and
// returning the new values for all N particles. New locations have to be shifted right after
// they are computed, because the new forces have to be calculated with all particles in the box
std::tuple<std::array<Vector2d, N>, std::array<Vector2d, N>, std::array<Vector2d, N>>
Verlet_step(const std::array<Vector2d, N> &loc_0, const std::array<Vector2d, N> &vel_0,
            const std::array<Vector2d, N> &F_0, double t_0, double h) {
    std::array<Vector2d, N> loc_1, vel_1;

    for (size_t i = 0; i < N; ++i) {
        loc_1[i] = (loc_0[i] + h * vel_0[i] + h * h * F_0[i] / (2 * M))
                       .unaryExpr([](double xi) { return xi - L * floor(xi / L); });
    }

    std::array<Vector2d, N> F_1;
    F_1 = calcForces(loc_1, t_0 + h);

    for (size_t i = 0; i < N; ++i) {
        vel_1[i] = vel_0[i] + h * (F_0[i] + F_1[i]) / (2 * M);
    }

    return std::make_tuple(loc_1, vel_1, F_1);
}

// lattice-like initialization with random velocities
// center of mass velocity is set to 0, temperature is scaled to match given T_0
std::tuple<std::array<Vector2d, N>, std::array<Vector2d, N>> initialize(double T_0) {
    int ppa = static_cast<int>(sqrt(N));
    double dist = L / (2 * ppa);

    if (pow(ppa, 2) != N) {
        std::cerr << "initialize: sqrt(N) should be an integer! EXIT" << std::endl;
        exit(EXIT_FAILURE);
    }

    namespace rand = boost::random;

    rand::mt11213b generator;                               // Mersenne Twister Generator
    rand::uniform_real_distribution<> normdist(-1.0, 1.0);  // Distribution
    rand::variate_generator<rand::mt11213b &, rand::uniform_real_distribution<>> norm_rnd(generator, normdist);

    std::array<Vector2d, N> loc, vel;
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
    double scale = sqrt(T_0 / temp);
    for (Vector2d &v : vel) v *= scale;

    Vector2d cm(CM_vel(vel));
    std::cout << "Init temperature: " << Temp(vel) << "\nInit CM velocity: ";
    std::cout << cm.norm() << std::endl;

    return std::make_tuple(loc, vel);
}

// function to get a traceback from python in case of an error
std::string extractPythonException() {
    namespace py = boost::python;

    PyObject *exc, *val, *tb;
    PyErr_Fetch(&exc, &val, &tb);
    PyErr_NormalizeException(&exc, &val, &tb);
    py::handle<> hexc(exc), hval(py::allow_null(val)), htb(py::allow_null(tb));

    if (!hval) {
        return py::extract<std::string>(py::str(hexc));
    } else {
        py::object traceback(py::import("traceback"));
        py::object format_exception(traceback.attr("format_exception"));
        py::object formatted_list(format_exception(hexc, hval, htb));
        py::object formatted(py::str("").join(formatted_list));
        return py::extract<std::string>(formatted);
    }
}

// sampling function to add all particles to the location and velocity samples
template <typename T1, typename T2>
void add_to_sample(const std::array<Eigen::Matrix<T1, 2, 1>, N> &vec, boost::multi_array<T2, 3> &sample, size_t at) {
    size_t i = 0;
    for (const Eigen::Matrix<T1, 2, 1> &r : vec) {
        sample[at][i][0] = static_cast<T2>(r[0]);
        sample[at][i][1] = static_cast<T2>(r[1]);
        i++;
    }
}

int main(int argc, char const *argv[]) {
    if (argc != 2) {
        std::cerr << "Function only takes one Argument! Please only pass a Value for T_0! EXIT" << std::endl;
        exit(EXIT_FAILURE);
    }
    double start = omp_get_wtime();

    // setting constants
    const double h = 0.01;

    const double t_0 = 0;
    const double t_max = 1e4;

    const double T_0 = boost::lexical_cast<double>(argv[1]);
    // const double T_0 = 0.0;
    // Init temperature varies from mean temperature after equilibration
    // const double T_0 = 1.000;  // yiedls T_mean = 1.4467
    // const double T_0 = 0.800;  // yiedls T_mean = 1.2614
    // const double T_0 = 0.600;  // yiedls T_mean = 1.0957
    // const double T_0 = 0.500;  // yiedls T_mean = 1.0072
    // const double T_0 = 0.300;  // yiedls T_mean = 0.8508
    // const double T_0 = 0.250;  // yiedls T_mean = 0.8183
    // const double T_0 = 0.001;  // yiedls T_mean = 0.6588
    // const double T_0 = 1e-05;  // yiedls T_mean = 0.6519
    // const double T_0 = 1e-10;  // yiedls T_mean = 0.6474

    // steps and number of samples
    const size_t samples = 10000;

    const size_t steps = static_cast<size_t>((t_max - t_0) / h + 0.5);
    const size_t sbs = static_cast<size_t>(steps * 1. / samples);  // steps between samples

    if (samples > steps) {
        std::cerr << "Number of samples should be less or equal number of steps! " << std::endl;
        exit(EXIT_FAILURE);
    }

    // arrays for storage of sample data
    using dtype = float;
    boost::multi_array<dtype, 3> loc_sample(boost::extents[samples][N][2]);
    boost::multi_array<dtype, 3> vel_sample(boost::extents[samples][N][2]);
    std::array<dtype, samples> E_kin_sample;
    std::array<dtype, samples> E_pot_sample;
    std::array<dtype, samples> cm_vel_sample;
    std::array<float, bins> radial_sample;

    // initialization
    std::array<std::array<Vector2d, N>, 2> locations, velocities, forces;
    std::tie(locations[0], velocities[0]) = initialize(T_0);
    forces[0] = calcForces(locations[0], t_0);
    // for (const Vector2d &f : forces[0]) std::cout << f << std::endl;

    // sample init status
    add_to_sample(locations[0], loc_sample, 0);
    add_to_sample(velocities[0], vel_sample, 0);

    E_kin_sample[0] = static_cast<dtype>(N * Temp(velocities[0]));
    cm_vel_sample[0] = static_cast<dtype>(CM_vel(velocities[0]).norm());
    E_pot_sample[0] = static_cast<dtype>(V_LJ(locations[0], t_0));

    // start verlet algorithm
    std::cout << "Verlet in progress:  " << std::flush;
    size_t at = 1;

    for (size_t t = 1; t < steps; t++) {
        std::tie(locations[t % 2], velocities[t % 2], forces[t % 2]) =
            Verlet_step(locations[(t + 1) % 2], velocities[(t + 1) % 2], forces[(t + 1) % 2], t_0 + t * h, h);

        if (t > t_aqui) update_radial(locations[t % 2], radial_sample, bins);

        if (t % sbs == 0) {
            add_to_sample(locations[t % 2], loc_sample, at);
            add_to_sample(velocities[t % 2], vel_sample, at);

            E_kin_sample[at]  = static_cast<dtype>(N * Temp(velocities[t % 2]));
            cm_vel_sample[at] = static_cast<dtype>(CM_vel(velocities[t % 2]).norm());
            E_pot_sample[at]  = static_cast<dtype>(V_LJ(locations[0], t_0));

            printf((10 * t < steps) ? "%2.1f%%\b\b\b\b" : "%2.1f%%\b\b\b\b\b", static_cast<double>(100 * t / steps));
            at++;
        }
    }

    // renormalization of g(r)
    for (size_t i = 0; i < bins; i++) {
        double dV = M_PI * (pow(i + 1, 2.0) - pow(i, 2.0)) * pow(cutoff / bins, 2.0);
        radial_sample[i] *= static_cast<float>(L * L / (dV * N * N * steps));
    }

    printf("%3.1f%%\ndone\n", 100.);
    std::cout << "seconds passed: " << omp_get_wtime() - start << std::endl;

    //––– plotting with python ––––––––––––––––––––––––––––––––––––––––––––––––
    namespace py = boost::python;
    namespace np = boost::numpy;

    try {
        Py_Initialize();
        np::initialize();

        // Retrieve the main module's namespace
        py::object global(py::import("__main__").attr("__dict__"));

        // Import neccessary modules
        py::exec("from __future__ import division \n"
                 "print 'Hello from Python!' \n"
                 "import numpy as np \n"
                 "print 'importing matplotlib...' \n"
                 "from matplotlib import pyplot as plt \n"
                 "from matplotlib import animation\n"
                 "plt.switch_backend('AGG') \n"
                 "plt.style.use('ggplot') \n"
                 "print 'done' \n"
                 "print 'plotting...' \n"
                 , global, global);

        global["file_ext"] = ".png";

        // Import variables and vectors
        global["L"] = L;
        global["N"] = N;
        global["T_0"] = T_0;
        global["t_0"] = t_0;
        global["t_max"] = t_max;
        global["t_aqui"] = t_aqui;
        global["sbs"] = sbs;
        global["samples"] = samples;
        global["cutoff"] = cutoff;
        global["bins"] = bins;

        // global["velocity_scale"] = 12;
        // global["force_scale"] = 1e-3;

        py::tuple strides_1 =
            static_cast<py::tuple>(py::eval("np.ndarray((samples, N, 2), np.float32).strides", global, global));

        global["loc_sample"] = np::from_data(loc_sample.data(), np::dtype::get_builtin<dtype>(),
                                             py::make_tuple(samples, N, 2), strides_1, py::object());
        global["vel_sample"] = np::from_data(vel_sample.data(), np::dtype::get_builtin<dtype>(),
                                             py::make_tuple(samples, N, 2), strides_1, py::object());
        global["E_kin_sample"] = np::from_data(E_kin_sample.data(), np::dtype::get_builtin<dtype>(),
                                               py::make_tuple(samples), py::make_tuple(sizeof(dtype)), py::object());
        global["E_pot_sample"] = np::from_data(E_pot_sample.data(), np::dtype::get_builtin<dtype>(),
                                               py::make_tuple(samples), py::make_tuple(sizeof(dtype)), py::object());
        global["cm_vel_sample"] = np::from_data(cm_vel_sample.data(), np::dtype::get_builtin<dtype>(),
                                                py::make_tuple(samples), py::make_tuple(sizeof(dtype)), py::object());
        global["radial_sample"] = np::from_data(radial_sample.data(), np::dtype::get_builtin<dtype>(),
                                                py::make_tuple(bins), py::make_tuple(sizeof(dtype)), py::object());

        // Launch some function in Python
        py::exec_file("plots.py", global, global);
    } catch (const py::error_already_set &) {
        std::cout << extractPythonException() << std::endl;
        exit(EXIT_FAILURE);
    }

    return 0;
}
