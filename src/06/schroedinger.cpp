#include <Eigen/Sparse>
#include <Eigen/LU>

#include <boost/python.hpp>
#include <boost/numpy.hpp>

#include <Python.h>

#include <iostream>
#include <functional>
#include <array>
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

int main() {
    const double L(10);  // computation window is [-L, +L] in both directions
    const double d_l(2);
    const double d_t(0.02);

    const int N(static_cast<int>(2 * L / d_l + 0.5));

    Eigen::SparseMatrix<std::complex<double>> H(N, N);
    Eigen::MatrixXcd S_p(N, N), S_n(N, N);

    H.reserve(3 * N - 2);

    for (int i = 0; i < N; i++) {
        H.insert(i, i) = 2. / (d_l * d_l) + d_l * d_l * i * i;

        S_p(i, i) = (1 + H.coeffRef(i, i).imag(), 0.5 * d_t * H.coeffRef(i, i).real());
        S_n(i, i) = (1 + H.coeffRef(i, i).imag(), -0.5 * d_t * H.coeffRef(i, i).real());
    }
    for (int i = 0; i < N - 1; i++) {
        H.insert(i + 1, i) = -1. / (d_l * d_l);
        H.insert(i, i + 1) = -1. / (d_l * d_l);

        S_p(i + 1, i) = (1 + H.coeffRef(i + 1, i).imag(), 0.5 * d_t * H.coeffRef(i + 1, i).real());
        S_p(i, i + 1) = (1 + H.coeffRef(i, i + 1).imag(), 0.5 * d_t * H.coeffRef(i, i + 1).real());
        S_n(i + 1, i) = (1 + H.coeffRef(i + 1, i).imag(), -0.5 * d_t * H.coeffRef(i + 1, i).real());
        S_n(i, i + 1) = (1 + H.coeffRef(i, i + 1).imag(), -0.5 * d_t * H.coeffRef(i, i + 1).real());
    }

    Eigen::MatrixXcd S = S_p.inverse() * S_n;
    std::cout << S << std::endl;

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
                 "import numpy as np \n",
                 global, global);

        // Import variables and vectors
        // global["L_x"] = L_x;
        // global["L_y"] = L_y;
        // global["N_x"] = N_x;
        // global["N_y"] = N_y;

        // py::tuple strides_A =
        // static_cast<py::tuple>(py::eval("np.ndarray((N_x, N_y), np.float64).strides", global, global));
        // np::dtype dtype = np::dtype::get_builtin<double>();

        // global["A"] = np::from_data(A.data(), dtype, py::make_tuple(N_x, N_y), strides_A, py::object());

        // Launch some function in Python
        // py::exec_file("plots.py", global, global);
    } catch (const py::error_already_set &) {
        std::cout << extractPythonException() << std::endl;
        exit(EXIT_FAILURE);
    }

    return 0;
}
