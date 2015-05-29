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



template<size_t N, size_t numcol>
void iterx(std::array<float,N> & arr, size_t transient, size_t steady, size_t colindex, std::function<double(double)> f){
	std::srand(static_cast<unsigned int>(N + colindex));
	double x1 = 0.5;
	double x2 = (static_cast<double>(rand())/static_cast<double>(RAND_MAX));
	for(size_t i = 0; i <transient+steady; i++){
		x1 = x2;
		x2 = f(x1);
		if (i < transient) continue;
		arr[(i - transient) * numcol + colindex] = static_cast<float>(x2) ;
	}
}

int main() {
  constexpr size_t steps = 500;
  constexpr size_t steady = 300;
  constexpr size_t N = steady * steps;

  const double rmin1 = 2.6;
  const double rmax1 = 4.0;
  const double dr1 = (rmax1 - rmin1) / steps;
  const double rmin2 = 1.6;
  const double rmax2 = 2.5;
  const double dr2 = (rmax2 - rmin2) / steps;
  const size_t transient = 100;

  std::array<float, N> data1, data2;

  double r1 = rmin1;
  double r2 = rmin2;
  for (size_t j = 0; j < steps; j++){
    r1 += dr1;
	r2 += dr2;
	 // good range for r is [2.6,4]
	iterx<N,steps>(data1, transient, steady, j, [r1](double x){return r1 * x * (1 - x); } );
	 // good range for r is [1.6,2.5]
	iterx<N,steps>(data2, transient, steady, j, [r2](double x){return r2 * x - pow(x, 3); } );
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
    py::exec(
        "from __future__ import division \n"
        "print 'Hello from Python!' \n"
        "import numpy as np \n",
        global, global);

    // Import variables and vectors
  global["steady"] = steady;
  global["rmin1"] = rmin1;
  global["dr1"] = dr1;
  global["rmin2"] = rmin2;
  global["dr2"] = dr2;
  global["steps"] = steps;

global["data1"] = np::from_data(
    data1.data(), np::dtype::get_builtin<float>(), py::make_tuple(N),
    py::make_tuple(sizeof(float)), py::object());

global["data2"] = np::from_data(
    data2.data(), np::dtype::get_builtin<float>(), py::make_tuple(N),
    py::make_tuple(sizeof(float)), py::object());

    // Launch some function in Python
    py::exec_file("plots.py", global, global);
  } catch (const py::error_already_set &) {
    std::cout << extractPythonException() << std::endl;
        exit(EXIT_FAILURE);
    }

    return 0;
}
