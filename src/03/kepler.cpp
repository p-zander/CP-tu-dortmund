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

  vector<VectorXd> r;
  vector<VectorXd> v;

  std::tie(r, v) =
      RK4_newton(r0, v0, 0, 100, 1000,
                 [](Vector3d r, double t) { return -r / pow(r.norm(), 3); });

  cout << "# x \t y \t z \t v_x \t v_y \t v_z \t Ldot \t Lx \t Ly \t Lz" << endl;

  for (auto Ir(r.begin()), Iv(v.begin()); Ir != r.end() && Iv != v.end();
       Ir++, Iv++) {
    Vector3d R(*Ir);
    Vector3d V(*Iv);
    Vector3d L(R.cross(V));
    Vector3d Ldot(R.cross(-R / pow(R.norm(), 3)));
    cout << R[0] << "\t" << R[1] << "\t" << R[2] << "\t" << V[0] << "\t" << V[1]
         << "\t" << V[2] << "\t" << Ldot.norm() << "\t" << L[0] << "\t" << L[1]
         << "\t" << L[2] << "\t" << endl;
  }
  return 0;
}