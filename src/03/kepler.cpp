#include "rungekutta.h"
#include <iostream>
#include <Eigen/Dense>
#include <functional>
#include <utility>
#include <vector>
#include <cmath>

using std::vector;
using std::cout;
using std::endl;
using Eigen::Vector3d;

int main() {
  Vector3d r0(1, 0, 0);
  Vector3d v0(0, 1.3, 0);

  double t0 = 0;
  double tN = 100;
  size_t steps = 10000;

  vector<VectorXd> r;
  vector<VectorXd> v;
  vector<VectorXd> ru;
  vector<VectorXd> vu;

  std::tie(r, v) =
      RK4_newton(r0, v0, t0, tN, steps,
                 [](Vector3d r, double t) { return -r / pow(r.norm(), 3); });
  std::tie(ru, vu) =
      RK4_newton(r.back(), -v.back(), t0, tN, steps,
                 [](Vector3d r, double t) { return -r / pow(r.norm(), 3); });

  cout << "# t \t x \t y \t z \t v_x \t v_y \t v_z \t Ldot \t Lx \t Ly \t Lz "
          "\t LRx "
          "\t LRy \t LRz \t Lx \t Ly \t Lz \t LRx \t LRy \t LRz \t ru_x \t "
          "ru_y \t ru_z \t vu_x \t vu_y \t vu_z" << endl;

  for (unsigned int i = 0; i <= steps; i++) {
    Vector3d R(r[i]);
    Vector3d V(v[i]);
    Vector3d L(R.cross(V));
    Vector3d Ldot(R.cross(-R / pow(R.norm(), 3)));
    Vector3d LR(V.cross(L) - R / R.norm());
    cout << i*(tN - t0) / steps << "\t" << R[0] << "\t" << R[1] << "\t" << R[2]
         << "\t" << V[0] << "\t" << V[1] << "\t" << V[2] << "\t" << Ldot.norm()
         << "\t" << L[0] << "\t" << L[1] << "\t" << L[2] << "\t" << LR[0]
         << "\t" << LR[1] << "\t" << LR[2] << "\t" << ru.back()[0] << "\t"
         << ru.back()[1] << "\t" << ru.back()[2] << "\t" << vu.back()[0] << "\t"
         << vu.back()[1] << "\t" << vu.back()[2] << endl;
  }
  return 0;
}