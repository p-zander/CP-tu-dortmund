#include <Eigen/Core>
#include <functional>
#include <vector>
#include <cmath>
#include <iostream>

#include "rungekutta.h"

using std::vector; using std::cout; using std::endl; using Eigen::Vector3d;

int main() {
    Vector3d r0(1,0,0);
    Vector3d v0(0,1,0);

    vector<VectorXd> r;
    vector<VectorXd> v;

    size_t steps = 100;
    double t0 = 0;
    double tN = 10;
    double h  = (tN-t0)/steps;

    std::tie(r, v) = RK4_newton(r0, v0, t0, tN, steps, [](Vector3d r, double t){ return -r; });

    cout << "# t \t x \t y \t z \t v_x \t v_y \t v_z" << endl;

    for (size_t i = 0; i <= steps; i++) { 
        cout <<   i*h   << "\t"
             << r[i][0] << "\t" << r[i][1] << "\t" << r[i][2] << "\t"
             << v[i][0] << "\t" << v[i][1] << "\t" << v[i][2] << "\t" << endl;
    }

    return 0;
}