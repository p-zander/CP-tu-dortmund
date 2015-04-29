#include <Eigen/Core>
#include <functional>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

#include "./rungekutta.h"

using std::vector;
using std::cout;
using std::endl;
using std::ofstream;
using Eigen::Vector3d;

int main() {
    Vector3d r0(1, 0, 0);
    double t0 = 0;
    double tN = 8*3.14159265358;
    
    Vector3d v0(0, 0, 0);
    size_t steps = 24;  // h = 1.0472
    double h = (tN - t0) / steps;

    vector<VectorXd> r, v;

    std::tie(r, v) = RK4_newton(r0, v0, t0, tN, steps, [](Vector3d r, double) { return -r; });

    ofstream file1("2a.txt");
    file1 << "# t \t x \t y \t z \t v_x \t v_y \t v_z" << endl;
    for (size_t i = 0; i <= steps; i++) {
        file1 <<   i*h   << "\t"
              << r[i][0] << "\t" << r[i][1] << "\t" << r[i][2] << "\t"
              << v[i][0] << "\t" << v[i][1] << "\t" << v[i][2] << "\t" << endl;
    }

    v0 = Vector3d(0, 1, 0);
    steps = 88;  // h = 0.2856
    h = (tN - t0) / steps;

    r.clear();
    v.clear();
    std::tie(r, v) = RK4_newton(r0, v0, t0, tN, steps, [](Vector3d r, double) { return -r; });

    ofstream file2("2b.txt");

    file2 << "# t \t x \t y \t z \t v_x \t v_y \t v_z" << endl;

    for (size_t i = 0; i <= steps; i++) {
        file2 <<   i*h   << "\t"
              << r[i][0] << "\t" << r[i][1] << "\t" << r[i][2] << "\t"
              << v[i][0] << "\t" << v[i][1] << "\t" << v[i][2] << "\t" << endl;
    }

    return 0;
}
