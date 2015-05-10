#include <boost/python.hpp>
#include <iostream>
#include <Python.h>
#include <utility>
#include <vector>

#include "./utils.h"

namespace py = boost::python;
using std::cout; using std::endl; using std::vector; using std::make_pair;

int main(int argc, char const *argv[]) {
    vector<int> x {1, 2, 5, 6, 10};
    vector<int> y {10, 5, 6, 3, -1};

    try {
        Py_Initialize();

        // Retrieve the main module.
        py::object main = py::import("__main__");

        // Retrieve the main module's namespace
        py::object global(main.attr("__dict__"));

        utils::to_numpy(make_pair(x, y), make_pair("x", "y"), &global, &global);

        // Launch some function in Python.
        py::exec("print 'Hello from Python!' \n"
                 "from matplotlib import pyplot as plt \n"
                 "import numpy as np \n"
                 "plt.plot(x, y) \n"
                 "plt.savefig('fig.pdf') \n"
                 "number = 42 \n"
                 , global, global);

        cout << py::extract<int>(global["number"]) << " = 42" << endl;

        Py_Finalize();
    }catch(const py::error_already_set&) {
        cout << pyutils::extractPythonException() << endl;
        exit(EXIT_FAILURE);
    }

    exit(EXIT_SUCCESS);
}
