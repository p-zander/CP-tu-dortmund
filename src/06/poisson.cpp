#include <Eigen/Core>

#include <boost/python.hpp>
#include <boost/numpy.hpp>

#include <Python.h>

#include <iostream>
#include <string>

using Eigen::Matrix;
using Eigen::MatrixXd;

template <typename T> void jacobi_iter_step(Matrix<T, -1, -1> &A, const Matrix<T, -1, -1> &rho, const double d_l) {
    if (A.cols() != rho.cols() || A.rows() != rho.rows()) {
        std::cerr << "Matrix 'A' and matrix 'rho': Dimension mismatch. EXIT" << std::endl;
        exit(EXIT_FAILURE);
    }

    for (int i = 1; i < A.cols() - 1; i++) {
        for (int j = 1; j < A.rows() - 1; j++) {
            A(j, i) = 1. / 4 * (A(j - 1, i) + A(j + 1, i) + A(j, i - 1) + A(j, i + 1) + d_l * d_l * rho(j, i));
        }
    }
}

// faulty function
// template <typename T>
// void gauss_seidel_iter_step(Matrix<T, -1, -1> &A, const Matrix<T, -1, -1> &rho, const double d_l) {
//     if (A.cols() != rho.cols() || A.rows() != rho.rows()) {
//         std::cerr << "Matrix 'A' and matrix 'rho': Dimension mismatch. EXIT" << std::endl;
//         exit(EXIT_FAILURE);
//     }

//     for (int i = 1; i < A.cols() - 2; i += 2) {
//         for (int j = 1; j < A.rows() - 2; j += 2) {
//             A(j, i) = 1. / 4 * (A(j - 1, i) + A(j + 1, i) + A(j, i - 1) + A.row(j).col(i + 1) + d_l * d_l * rho(j,
//             i));
//         }
//     }
//     for (int i = 2; i < A.cols() - 1; i += 2) {
//         for (int j = 2; j < A.rows() - 1; j += 2) {
//             A(j, i) = 1. / 4 * (A(j - 1, i) + A(j + 1, i) + A(j, i - 1) + A(j, i + 1) + d_l * d_l * rho(j, i));
//         }
//     }
// }

template <typename T, class iter_step_F>
void do_some_iteration(Matrix<T, -1, -1> &A, const iter_step_F &iter_step, const int start = 10, int step = 5,
                       const double eps = 1e-5) {
    Matrix<T, -1, -1> copy_A(A.cols(), A.rows());
    Matrix<T, -1, -1> EPS(A.cols(), A.rows());

    EPS = eps * Matrix<T, -1, -1>::Ones(A.cols(), A.rows());

    for (int i = 0; i < start - 1; i++) iter_step(A);
    copy_A = A;
    iter_step(A);

    int N = start;
    while (((A - copy_A).cwiseAbs().array() > EPS.array()).any()) {
        N += step;

        for (int i = 0; i < step - 1; i++) iter_step(A);
        copy_A = A;
        iter_step(A);
    }

    std::cout << "Steps needed: " << N << std::endl;
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
    rho(rho.rows() / 2, rho.cols() / 2) = 1;
}

// suits the assignment only with even N
template <typename T> void init_e(Matrix<T, -1, -1> &__attribute__((unused)), Matrix<T, -1, -1> &rho, int N) {
    if (N % 2 != 0) {
        std::cerr << "init_e: integer N should be an even number. EXIT" << std::endl;
        exit(EXIT_FAILURE);
    }

    for (int i = 1; i <= N; i++) {
        for (int j = 1; j <= N; j++) {
            int charge = (i % 2) ? -1 : 1;
            charge *= (j % 2) ? -1 : 1;

            rho((rho.rows() * i / (N + 1)), (rho.cols() * j / (N + 1))) = charge;
        }
    }
}

int main(int argc, char const *argv[]) {
    const double L(1);
    const double d_l(0.05);

    const unsigned int N(static_cast<unsigned int>(L / d_l + 0.5));

    MatrixXd A(N, N);
    MatrixXd rho(N, N);

    std::string part(argv[1]);

    if (argc > 1 && part == "b") {
        init_b(A, rho);
    } else {
        if (argc > 1 && part == "c") {
            init_c(A, rho);
        } else {
            if (argc > 2 && part == "e") {
                init_e(A, rho, std::stoi(std::string(argv[2])));
            } else {
                std::cerr << "input: 'b', 'c' or 'e N' where N is an even integer. EXIT" << std::endl;
                exit(EXIT_FAILURE);
            }
        }
    }

    do_some_iteration(A, [&rho, d_l](MatrixXd &B) { jacobi_iter_step(B, rho, d_l); });

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
        global["part"] = part;
        global["L"] = L;
        global["N"] = N;

        py::tuple strides_A =
            static_cast<py::tuple>(py::eval("np.ndarray((N, N), np.float64).strides", global, global));
        np::dtype dtype = np::dtype::get_builtin<double>();

        global["A"] = np::from_data(A.data(), dtype, py::make_tuple(N, N), strides_A, py::object());

        // Launch some function in Python
        py::exec_file("2_plots.py", global, global);
    } catch (const py::error_already_set &) {
        std::cout << extractPythonException() << std::endl;
        exit(EXIT_FAILURE);
    }

    return 0;
}
