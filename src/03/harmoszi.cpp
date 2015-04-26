#include <Eigen/Core>
#include <functional>
#include <vector>
#include <cmath>
#include <iostream>

#include "rungekutta.h"

using std::vector; using std::cout; using std::endl; using Eigen::Vector3d;

int main() {
    Vector3d r0(1,0,0);
    Vector3d v0(0,0,0);

    vector<VectorXd> r = RK4_newton(r0, v0, 0, 10, 100,[](Vector3d r, double t){ return -r; });

    for(const auto& x : r) 
        cout << x[0] << "\t" << x[1] << "\t" << x[2] << "\t" << endl;

    return 0;
}