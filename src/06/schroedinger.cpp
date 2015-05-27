#include <Eigen/Sparse>
#include <Eigen/LU>

#include <boost/python.hpp>
#include <boost/numpy.hpp>
#include <boost/multi_array.hpp>

#include <Python.h>

#include <iostream>
#include <functional>
#include <array>
#include <string>

using cd = std::complex<double>;
using VecXcd = Eigen::VectorXcd;
using matXcd = Eigen::MatrixXcd;
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
  const double L{10};  // computation window is [-L, +L] in both directions
  const double d_l{.1};
  const double d_t{0.02};
  const double xi0{1};
  const double si{1};
  const size_t t{500};

  const int N(1 + static_cast<int>(2 * L / d_l + 0.5));

  matXcd H = matXcd::Identity(N, N) * 2;

  H.topRightCorner(N - 1, N - 1) += matXcd::Identity(N - 1, N - 1);
  H.bottomLeftCorner(N - 1, N - 1) += matXcd::Identity(N - 1, N - 1);
  H *= -1. / pow(d_l, 2);
  H += pow(d_l,2) * Eigen::VectorXd::LinSpaced(N, -L, L).cwiseAbs2().cast<cd>().asDiagonal();
  matXcd S_p = matXcd::Identity(N, N) + std::complex<double>(0, 1) * 0.5 * d_t * H;
  matXcd S_n = matXcd::Identity(N, N) - std::complex<double>(0, 1) * 0.5 * d_t * H;
  matXcd S = S_p.inverse() * S_n;

  VecXcd psi0 = Eigen::VectorXd::LinSpaced(N, -L, L)
                    .unaryExpr([xi0, si](double xi) {
                      return exp(-pow(xi - xi0, 2) / (4 * si));
                    })
                    .cast<cd>()
                    .normalized();

  std::vector<VecXcd> psi(t, VecXcd(N));
  psi[0] = psi0;

  for (std::vector<VecXcd>::iterator it = psi.begin() + 1; it != psi.end();++it) *it = S * *(it - 1);

  boost::multi_array<float, 2> psi2(boost::extents[N][t]);
  for (long i = 0; i < N; i++) {
    for (long j = 0; j < static_cast<long>(t); j++) {
      psi2[i][j] = static_cast<float>(norm(psi[static_cast<size_t>(j)][i]));
    }
  }

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
        global["N"] = N;
        global["t"] = t;
        // global["N_x"] = N_x;
        // global["N_y"] = N_y;

        py::tuple strides_1 = static_cast<py::tuple>(py::eval(
            "np.ndarray((N, t), np.float32).strides", global, global));

        global["psi2"] = np::from_data(
            psi2.data(), np::dtype::get_builtin<float>(),
            py::make_tuple(N, t), strides_1, py::object());

        // Launch some function in Python
        py::exec_file("plots2.py", global, global);
    } catch (const py::error_already_set &) {
        std::cout << extractPythonException() << std::endl;
        exit(EXIT_FAILURE);
    }

    return 0;
}
