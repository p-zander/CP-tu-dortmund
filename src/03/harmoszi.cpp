#include <Eigen/Core>
#include <utility>
#include <vector>
#include <cmath>
#include <fstream>

#include "./rungekutta.h"

using std::vector;
using std::endl;
using std::ofstream;
using Eigen::Vector3d;

int main() {
    Vector3d r0(1, 0, 0);
    double t0 = 0;
    double tN = 4*M_PI;

    Vector3d v0(0, 0, 0);
    size_t steps = 12;  // h = 1.0472
    double h = (tN - t0) / steps;

    vector<VectorXd> r, v;
    /* can't be vector<Vector3d> because casting vector<x> to vector<y> can't be
    done implicitly even if x can be implicitly casted to y */

    std::tie(r, v) = RK4_newton(r0, v0, t0, tN, steps, [](Vector3d r, double) { return -r; });

    ofstream file1("2a.txt");
    file1 << "# t \t x \t y \t z \t v_x \t v_y \t v_z" << endl;
    for (size_t i = 0; i <= steps; i++) {
        file1 <<   i*h   << "\t"
              << r[i][0] << "\t" << r[i][1] << "\t" << r[i][2] << "\t"
              << v[i][0] << "\t" << v[i][1] << "\t" << v[i][2] << "\t" << endl;
    }

    v0 = Vector3d(0, 0.5, 0);
    steps = 52;  // h = 0.2618
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
