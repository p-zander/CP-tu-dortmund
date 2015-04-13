#include <fstream>
#include <cmath>
#include <Eigen/Dense>
#include <iostream>
using namespace std;

// Energy E
static double E(int N, double theta, bool antiferro = false) {
  theta *= M_PI / 180;
  double sum = 0;
  Eigen::Vector2d R(0, 0);
  Eigen::Vector2d m(sin(theta), cos(theta));
  Eigen::Vector2d n(0, 1);

  for (int a = -N; a <= N; ++a) {
    for (int b = -N; b <= N; ++b) {
      if (!(a == 0 && b == 0)) {
        R[0] = a;
        R[1] = b;
        if (antiferro) n[1] = pow(-1, abs(a + b));
        sum += (-3 * R.dot(m) * R.dot(n)) / pow(R.norm(), 5) +
               m.dot(n) / pow(R.norm(), 3);
      }
    }
  }

  return 1e-6 * sum;
}

// T via derivation of E
static double T(int N, double theta, double h, bool antiferro = false) {
  return -(E(N, theta - 2 * h, antiferro) - 8 * E(N, theta - h, antiferro) +
           8 * E(N, theta + h, antiferro) - E(N, theta + 2 * h, antiferro)) /
         (12 * h);
}

// T via m cross product
static double T2(int N, double theta, bool antiferro = false) {
  theta *= M_PI / 180;
  Eigen::Vector3d R(0, 0, 0);
  Eigen::Vector3d m(sin(theta), cos(theta), 0);
  Eigen::Vector3d n(0, 1, 0);
  Eigen::Vector3d B(0, 0, 0);
  for (int a = -N; a <= N; ++a) {
    for (int b = -N; b <= N; ++b) {
      if (!(a == 0 && b == 0)) {
        R[0] = b;
        R[1] = a;
        if (antiferro) n[1] = pow(-1, abs(a + b));
        B += 3 * R.dot(n) / pow(R.norm(), 5) * R - n / pow(R.norm(), 3);
      }
    }
  }
  return 1e-6 * (m.cross(B))[2];
}

// simple function to test correctness of T and T2
static double T3(int N, double theta, bool antiferro = false) {
  theta *= M_PI / 180;
  double sum = 0;
  double r = 0;
  short anti =1;
  for (int a = -N; a <= N; ++a) {
    for (int b = -N; b <= N; ++b) {
      if (!(a == 0 && b == 0)) {
        r = sqrt(a * a + b * b);
        if (antiferro) anti=pow(-1,abs(a+b));
        sum +=
            anti*(3 * b * (a * cos(theta) - b * sin(theta)) + r * r * sin(theta)) /
            pow(r, 5);
      }
    }
  }
  return 1e-6*sum;
}

int main() {
  for (int i = 0; i <= 360; ++i) {
    cout << i << "  " << T(2, i, 5) << "  " << T2(2, i) << " " << T3(2, i)
         << endl;
  }
  // ofstream file;
  // file.open("2NT_ferro.txt");
  // for (int i = 0; i <= 360; ++i) {
  //   file << i << " " << E(2, i, true) << " " << T(2, i,1) << endl;
  // }
  // file.close();
  return 0;
}
