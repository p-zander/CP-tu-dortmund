#include <Eigen/Core>

#include <boost/python.hpp>
#include <boost/numpy.hpp>

#include <Python.h>

#include <iostream>
#include <functional>
#include <tuple>
#include <string>

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

std::tuple<double, double> regula_falsi(const std::function<double(double)> &f, double x_0, double y_0) {
    double z_1 = (x_0 * f(y_0) - y_0 * f(x_0)) / (f(y_0) - f(x_0));
    if (f(x_0) * f(z_1) < 0)
        return std::make_tuple(x_0, z_1);
    else
        return std::make_tuple(z_1, y_0);
}

double regula_richti(const std::function<double(double)> &f, double x_0, double y_0, double eps = 1e-15) {
    while (y_0 - x_0 > eps) std::tie(x_0, y_0) = regula_falsi(f, x_0, y_0);
    return (y_0 + x_0) / 2;
}

std::tuple<double, double> bisection_step(const std::function<double(double)> &f, double x_0, double y_0) {
    double z_1 = (x_0 + y_0) / 2;
    if (f(x_0) * f(z_1) < 0)
        return std::make_tuple(x_0, z_1);
    else
        return std::make_tuple(z_1, y_0);
}

double bisection(const std::function<double(double)> &f, double x_0, double y_0, double eps = 1e-15) {
    while (y_0 - x_0 > eps) std::tie(x_0, y_0) = bisection_step(f, x_0, y_0);
    return (y_0 + x_0) / 2;
}

int main() {
    const double r_inf = 3.569946;
    unsigned int N = 1000;

    auto f = [](double r, double x) { return r * x * (1 - x); };

    std::function<double(double, double, unsigned int)> n_th_iter;
    n_th_iter = [&n_th_iter, f](double r, double x, unsigned int n) mutable -> double {
        if (n == 1)
            return f(r, x);
        else
            return f(r, n_th_iter(r, x, n - 1));
    };

    double x = 0.5;
    auto g_0 = [x, n_th_iter](double r) { return 0.5 - n_th_iter(r, x, 1); };
    auto g_1 = [x, n_th_iter](double r) { return 0.5 - n_th_iter(r, x, 2); };
    auto g_2 = [x, n_th_iter](double r) { return 0.5 - n_th_iter(r, x, 4); };
    auto g_3 = [x, n_th_iter](double r) { return 0.5 - n_th_iter(r, x, 8); };

    Eigen::VectorXd r_vec = Eigen::VectorXd::LinSpaced(N, 0, r_inf);
    Eigen::VectorXd g_0_vec = r_vec.unaryExpr(g_0);
    Eigen::VectorXd g_1_vec = r_vec.unaryExpr(g_1);
    Eigen::VectorXd g_2_vec = r_vec.unaryExpr(g_2);
    Eigen::VectorXd g_3_vec = r_vec.unaryExpr(g_3);

    double R_0 = bisection(g_0, 0, 3);
    double R_1 = bisection(g_1, 3, 3.3);
    double R_2 = bisection(g_3, 3.3, 3.5);
    double R_3 = bisection(g_3, 3.5, r_inf);

    std::cout.precision(12);
    std::cout << "Nullstellen:\n" << R_0 << "  " << R_1 << "  " << R_2 << "  " << R_3 << std::endl;
    std::cout << "Feigenbaum:  " << (R_2 - R_1) / (R_3 - R_2) << std::endl;

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
                 "from matplotlib import pyplot as plt \n"
                 "plt.style.use('ggplot') \n",
                 global, global);

        // Import variables and vectors
        global["r_inf"] = r_inf;

        global["r"] = np::from_data(r_vec.data(), np::dtype::get_builtin<double>(), py::make_tuple(N),
                                    py::make_tuple(sizeof(double)), py::object());
        global["g_0"] = np::from_data(g_0_vec.data(), np::dtype::get_builtin<double>(), py::make_tuple(N),
                                      py::make_tuple(sizeof(double)), py::object());
        global["g_1"] = np::from_data(g_1_vec.data(), np::dtype::get_builtin<double>(), py::make_tuple(N),
                                      py::make_tuple(sizeof(double)), py::object());
        global["g_2"] = np::from_data(g_2_vec.data(), np::dtype::get_builtin<double>(), py::make_tuple(N),
                                      py::make_tuple(sizeof(double)), py::object());
        global["g_3"] = np::from_data(g_3_vec.data(), np::dtype::get_builtin<double>(), py::make_tuple(N),
                                      py::make_tuple(sizeof(double)), py::object());

        // Launch some function in Python
        py::exec("plt.axhline(0, color='k') \n"
                 "plt.plot(r, g_0, label=r'$g_0$') \n"
                 "plt.plot(r, g_1, label=r'$g_1$') \n"
                 "plt.plot(r, g_2, label=r'$g_2$') \n"
                 "plt.plot(r, g_3, label=r'$g_3$') \n"
                 "plt.xlim(0, 4) \n"
                 "plt.ylim(-0.4, 0.6) \n"
                 "plt.legend(loc='best') \n"
                 "plt.savefig('2_a.pdf') \n",
                 global, global);
    } catch (const py::error_already_set &) {
        std::cout << extractPythonException() << std::endl;
        exit(EXIT_FAILURE);
    }

    return 0;
}
