#include <boost/python.hpp>
#include <iostream>
#include <Python.h>
#include <utility>
#include <vector>

#include "./utils.h"

using boost::python::object; using boost::python::exec; using boost::python::extract;
using std::cout; using std::endl; using std::vector; using std::make_pair;

int main(int argc, char const *argv[]) {
    vector<int> x {1, 2, 5, 6, 10};
    vector<int> y {10, 5, 6, 3, -1};

    try {
        Py_Initialize();

        // Retrieve the main module.
        object main = boost::python::import("__main__");

        // Retrieve the main module's namespace
        object global(main.attr("__dict__"));

        pyutils::to_numpy(make_pair(x, y), make_pair("x", "y"), &global, &global);

        // Launch some function in Python.
        exec("print 'Hello from Python!' \n"
             "from matplotlib import pyplot as plt \n"
             "import numpy as np \n"
             // "y = np.array([]) \n"
             // "print np.append(y,[11]) \n"
             "plt.plot(x, y) \n"
             "plt.savefig('fig.pdf') \n"
             "number = 42 \n"
             , global, global);

        cout << extract<int>(global["number"]) << " = 42" << endl;
    }catch(const boost::python::error_already_set& e) {
        cout << "Python Error. EXIT" << endl;  // No python Traceback though
        exit(EXIT_FAILURE);
    }

    return 0;
}
