#include <boost/python.hpp>
#include <iostream>
#include <Python.h>
#include <utility>
#include <vector>
#include <string>
#include <sstream>

#ifndef PYUTILS_H
#define PYUTILS_H

namespace utils {
template <typename T1, typename T2>
std::vector<std::pair<T1, T2>> to_vector_of_pair(const std::pair<std::vector<T1>, std::vector<T2>> &p) {
    assert(p.first.size() == p.second.size());

    std::vector<std::pair<T1, T2>> v;

    for (long i = 0; i < p.first.size(); i++) v.push_back(std::make_pair(p.first[i], p.second[i]));

    return v;
}

template <typename T1, typename T2>
std::pair<std::vector<T1>, std::vector<T2>> to_pair_of_vectors(const std::vector<std::pair<T1, T2>> &p) {
    std::vector<T1> v1;
    std::vector<T2> v2;

    for (const std::pair<T1, T2> &vals : p) {
        v1.push_back(vals.first);
        v2.push_back(vals.second);
    }

    return std::make_pair(v1, v2);
}

template <typename T1, typename T2> void to_numpy(const std::pair<std::vector<T1>, std::vector<T2>> &p,
                                                   const std::pair<std::string, std::string> &n,
                                                   boost::python::object *global, boost::python::object *local,
                                                   bool showcommands = false) {
    assert(p.first.size() == p.second.size());
    std::stringstream s;
    s << "import numpy as np" << std::endl << n.first << " = np.array([])" << std::endl;
    s << "import numpy as np" << std::endl << n.second << " = np.array([])" << std::endl;

    std::vector<std::pair<T1, T2>> v(to_vector_of_pair(p));

    s << std::scientific;
    for (const std::pair<T1, T2> &vals : v)
        s << n.first << " = np.append(" << n.first << ", " << vals.first << ")" << std::endl << n.second
          << " = np.append(" << n.second << ", " << vals.second << ")" << std::endl;

    if (showcommands) std::cout << s.str() << std::endl;
    boost::python::exec(s.str().c_str(), *global, *local);
}

template <typename T>
void to_numpy(const std::vector<T> &vec, const std::string &name, boost::python::object *global,
              boost::python::object *local, bool showcommands = false) {
    std::stringstream s;
    s << "import numpy as np" << std::endl << name << " = np.array([])" << std::endl;

    s << std::scientific;
    for (const T &val : vec) s << name << " = np.append(" << name << ", " << val << ")" << std::endl;

    if (showcommands) std::cout << s.str() << std::endl;
    boost::python::exec(s.str().c_str(), *global, *local);
}

template <typename T,
          int D>
void
to_numpy(const std::vector<Eigen::Matrix<T, D, 1>> &vec, const std::string &name, boost::python::object *global,
         boost::python::object *local, bool showcommands = false) {
    std::stringstream s;
    s << "import numpy as np" << std::endl << name << " = np.array([])" << std::endl;

    int i = 0;
    s << std::scientific;
    for (const Eigen::Matrix<T, D, 1> &val : vec) {
        s << name << " = np.append(" << name << ", [";
        for (long i = 0; i < D; i++) s << val[i] << ", ";
        s << "])" << std::endl;
        i++;
    }

    s << name << " = " << name << ".reshape((" << i << ", " << D << "))" << std::endl;

    if (showcommands) std::cout << s.str() << std::endl;
    boost::python::exec(s.str().c_str(), *global, *local);
}

template <typename T,
          int D, size_t N>
void
to_numpy(const std::array<Eigen::Matrix<T, D, 1>, N> &vec, const std::string &name, boost::python::object *global,
         boost::python::object *local, bool showcommands = false) {
    std::stringstream s;
    s << "import numpy as np" << std::endl << name << " = np.array([])" << std::endl;

    s << std::scientific;
    for (const Eigen::Matrix<T, D, 1> &val : vec) {
        s << name << " = np.append(" << name << ", [";
        for (long i = 0; i < D; i++) s << val[i] << ", ";
        s << "])" << std::endl;
    }

    s << name << " = " << name << ".reshape((" << N << ", " << D << "))" << std::endl;

    if (showcommands) std::cout << s.str() << std::endl;
    boost::python::exec(s.str().c_str(), *global, *local);
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
}  // namespace utils
#endif
