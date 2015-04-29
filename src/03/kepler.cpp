#include "rungekutta.h"
#include <Eigen/Dense>
#include <functional>
#include <utility>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>

using std::cout;
using std::vector;
using std::endl;
using Eigen::Vector3d;
using std::ofstream;
using std::stringstream;

void a(Vector3d v, double tN, size_t steps) {
  Vector3d r0(1, 0, 0);
  double t0 = 0;

  vector<VectorXd> R1;
  vector<VectorXd> V1;

  std::tie(R1, V1) =
      RK4_newton(r0, v, t0, tN, steps,
                 [](Vector3d r, double t) { return -r / pow(r.norm(), 3); });

  stringstream filename;
  filename << "a_tn_" << tN << "_N_" << steps << "_vy_" << v[1] << ".txt";
  ofstream file1;
  file1.open(filename.str().c_str());

  file1 << "# t \t x1 \t y1 \t z1 \t" << endl;
  for (unsigned int i = 0; i <= steps; i++) {
    file1 << i*(tN - t0) / steps << "\t" << R1[i][0] << "\t" << R1[i][1] << "\t"
          << R1[i][2] << endl;
  }
  file1.close();
}

void bcd() {
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

  ofstream file2;

  file2.open("bcd.txt");
  file2 << "# t \t x \t y \t z \t v_x \t v_y \t v_z \t Ldot \t Lx \t Ly \t Lz "
           "\t LRx \t LRy \t LRz \t ru_x \t ru_y \t ru_z \t vu_x \t vu_y \t "
           "vu_z" << endl;

  for (unsigned int i = 0; i <= steps; i++) {
    Vector3d R(r[i]);
    Vector3d V(v[i]);
    Vector3d L(R.cross(V));
    Vector3d Ldot(R.cross(-R / pow(R.norm(), 3)));
    Vector3d LR(V.cross(L) - R / R.norm());

    file2 << i*(tN - t0) / steps << "\t" << R[0] << "\t" << R[1] << "\t" << R[2]
          << "\t" << V[0] << "\t" << V[1] << "\t" << V[2] << "\t" << Ldot.norm()
          << "\t" << L[0] << "\t" << L[1] << "\t" << L[2] << "\t" << LR[0]
          << "\t" << LR[1] << "\t" << LR[2] << "\t" << ru[i][0] << "\t"
          << ru[i][1] << "\t" << ru[i][2] << "\t" << vu[i][0] << "\t"
          << vu[i][1] << "\t" << vu[i][2] << endl;
  }
  file2.close();
}

void e(double alpha) {
  Vector3d r0(1, 0, 0);
  Vector3d v0(0, 1.3, 0);

  double t0 = 0;
  double tN = 100;
  size_t steps = 1000;

  vector<VectorXd> r;
  vector<VectorXd> v;

  std::tie(r, v) =
      RK4_newton(r0, v0, t0, tN, steps, [alpha](Vector3d r, double t) {
        return -r * alpha / pow(r.norm(), alpha + 2);
      });

  ofstream file3;
  stringstream filename;
  filename << "e_alpha_" << alpha << ".txt";
  file3.open(filename.str().c_str());
  file3 << "# t \t x \t y \t z \t LRx \t LRy \t LRz" << endl;

  for (unsigned int i = 0; i <= steps; i++) {
    Vector3d R(r[i]);
    Vector3d V(v[i]);
    Vector3d L(R.cross(V));
    Vector3d LR(V.cross(L) - R / R.norm());

    file3 << i*(tN - t0) / steps << "\t" << R[0] << "\t" << R[1] << "\t" << R[2]
          << "\t" << LR[0] << "\t" << LR[1] << "\t" << LR[2] << endl;
  }
  file3.close();
}

int main() {
  Vector3d v0(0, 1.3, 0), v1(0, 0.88, 0);

  a(v0, 37, 100);
  a(v0, 37, 25);
  a(v1, 37, 100);
  bcd();
  e(0.9);
  e(1.1);

  return 0;
}