#include <Eigen/Core>

#include <boost/python.hpp>
#include <boost/numpy.hpp>

#include <Python.h>

#include <iostream>
#include <array>
#include <string>

using Eigen::Matrix;
using Eigen::MatrixXd;

template <typename T>
void jacobi_iter_step(Matrix<T, -1, -1> &A, const Matrix<T, -1, -1> &rho, double d_x, double d_y) {
    if (A.cols() != rho.cols() || A.rows() != rho.rows()) {
        std::cerr << "Matrix 'A' and matrix 'rho': Dimension mismatch. EXIT" << std::endl;
        exit(EXIT_FAILURE);
    }

    for (int i = 1; i < A.cols() - 1; i++) {
        for (int j = 1; j < A.rows() - 1; j++) {
            A.row(j).col(i) = 1. / 4 * (A.row(j - 1).col(i) + A.row(j + 1).col(i) + A.row(j).col(i - 1) +
                                        A.row(j).col(i + 1) + d_x * d_y * rho.row(j).col(i));
        }
    }
}

template <typename T>
void gauss_seidel_iter_step(Matrix<T, -1, -1> &A, const Matrix<T, -1, -1> &rho, double d_x, double d_y) {
    if (A.cols() != rho.cols() || A.rows() != rho.rows()) {
        std::cerr << "Matrix 'A' and matrix 'rho': Dimension mismatch. EXIT" << std::endl;
        exit(EXIT_FAILURE);
    }

    for (int i = 1; i < A.cols() - 2; i += 2) {
        for (int j = 1; j < A.rows() - 2; j += 2) {
            A.row(j).col(i) = 1. / 4 * (A.row(j - 1).col(i) + A.row(j + 1).col(i) + A.row(j).col(i - 1) +
                                        A.row(j).col(i + 1) + d_x * d_y * rho.row(j).col(i));
        }
    }
    for (int i = 2; i < A.cols() - 1; i += 2) {
        for (int j = 2; j < A.rows() - 1; j += 2) {
            A.row(j).col(i) = 1. / 4 * (A.row(j - 1).col(i) + A.row(j + 1).col(i) + A.row(j).col(i - 1) +
                                        A.row(j).col(i + 1) + d_x * d_y * rho.row(j).col(i));
        }
    }
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
template <typename T>
void init_a(Matrix<T, -1, -1> &__attribute__((unused)), Matrix<T, -1, -1> &__attribute__((unused))) {}

template <typename T> void init_b(Matrix<T, -1, -1> &A, Matrix<T, -1, -1> &__attribute__((unused))) {
    A.row(A.cols() - 1) = Matrix<T, -1, -1>::Ones(1, A.rows());
}

template <typename T> void init_c(Matrix<T, -1, -1> &__attribute__((unused)), Matrix<T, -1, -1> &rho) {
    rho.col(static_cast<int>(rho.rows() * 1. / 2)).row(static_cast<int>(rho.cols() * 1. / 3)) =
        Matrix<T, -1, -1>::Ones(1, 1);
}

template <typename T> void init_e_1(Matrix<T, -1, -1> &__attribute__((unused)), Matrix<T, -1, -1> &rho) {
    rho.col(static_cast<int>(rho.rows() * 1. / 3)).row(static_cast<int>(rho.cols() * 2. / 3)) +=
        Matrix<T, -1, -1>::Ones(1, 1);
    rho.col(static_cast<int>(rho.rows() * 2. / 3)).row(static_cast<int>(rho.cols() * 1. / 3)) -=
        Matrix<T, -1, -1>::Ones(1, 1);
}

template <typename T> void init_e_2(Matrix<T, -1, -1> &__attribute__((unused)), Matrix<T, -1, -1> &rho) {
    rho.col(static_cast<int>(rho.rows() * 1. / 3)).row(static_cast<int>(rho.cols() * 2. / 3)) +=
        Matrix<T, -1, -1>::Ones(1, 1);
    rho.col(static_cast<int>(rho.rows() * 2. / 3)).row(static_cast<int>(rho.cols() * 1. / 3)) +=
        Matrix<T, -1, -1>::Ones(1, 1);
    rho.col(static_cast<int>(rho.rows() * 2. / 3)).row(static_cast<int>(rho.cols() * 2. / 3)) -=
        Matrix<T, -1, -1>::Ones(1, 1);
    rho.col(static_cast<int>(rho.rows() * 1. / 3)).row(static_cast<int>(rho.cols() * 1. / 3)) -=
        Matrix<T, -1, -1>::Ones(1, 1);
}

int main() {
    const double L_x(1);
    const double L_y(1);
    const double d_x(0.05), d_y(0.05);

    const unsigned int N_x(static_cast<unsigned int>(L_x / d_x));
    const unsigned int N_y(static_cast<unsigned int>(L_y / d_y));

    MatrixXd A(N_y, N_x);
    MatrixXd rho(N_y, N_x);

    init_e_2(A, rho);

    for (int i = 0; i < 100; i++) {
        jacobi_iter_step(A, rho, d_x, d_y);
    }

    // std::cout << A << std::endl;

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
                 // "plt.switch_backend('AGG') \n"
                 "plt.style.use('ggplot') \n"
                 "print 'done' \n"
                 "print 'plotting...' \n",
                 global, global);

        global["file_ext"] = ".pdf";

        // Import variables and vectors
        global["L_x"] = L_x;
        global["L_y"] = L_y;
        global["N_x"] = N_x;
        global["N_y"] = N_y;

        py::tuple strides_A =
            static_cast<py::tuple>(py::eval("np.ndarray((N_y, N_x), np.float64).strides", global, global));

        global["A"] = np::from_data(A.data(), np::dtype::get_builtin<double>(), py::make_tuple(N_y, N_x), strides_A,
                                    py::object());
        // global["b"] = np::from_data(E_kin_sample.data(), np::dtype::get_builtin<dtype>(),
        // py::make_tuple(samples), py::make_tuple(sizeof(dtype)), py::object());

        py::exec("plt.imshow(A) \nplt.savefig('A.pdf')", global, global);

        // Launch some function in Python
        // py::exec_file("plots.py", global, global);
    } catch (const py::error_already_set &) {
        std::cout << extractPythonException() << std::endl;
        exit(EXIT_FAILURE);
    }

    return 0;
}
