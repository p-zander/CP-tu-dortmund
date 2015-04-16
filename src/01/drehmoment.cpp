#include <fstream>
#include <cmath>
#include <Eigen/Dense>
#include <sstream>
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
         (12 * h * M_PI / 180);
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
        R[0] = a;
        R[1] = b;
        if (antiferro) n[1] = pow(-1, abs(a + b));
        B += 3 * R.dot(n) / pow(R.norm(), 5) * R - n / pow(R.norm(), 3);
      }
    }
  }
  return -1e-6 * (m.cross(B))[2];
}

// Quick and dirty way to produce the needed text files.
int main() {
  int N[3] = {2, 5, 10};
  stringstream filename1;
  stringstream filename2;
  for (int a = 0; a < 2; ++a) {
    for (int b = 0; b < 3; ++b) {
      filename1.str("");
      filename2.str("");
      (a) ? filename1 << "E_Antiferro_" << N[b] << ".txt"
          : filename1 << "E_Ferro_" << N[b] << ".txt";
      (a) ? filename2 << "T_Antiferro_" << N[b] << ".txt"
          : filename2 << "T_Ferro_" << N[b] << ".txt";
      ofstream file1;
      ofstream file2;
      file1.open(filename1.str().c_str());
      file2.open(filename2.str().c_str());
      file1 << "theta  E" << endl;
      file2 << "theta  T(derived)  T(crosproduct)" << endl;
      for (int i = 0; i <= 360; ++i) {
        file1 << i << "  " << E(N[b], i, a) << endl;
        file2 << i << "  " << T(N[b], i, 1, a) << "  " << T2(N[b], i, a)
              << endl;
      }
      file1.close();
      file2.close();
    }
  }

  return 0;
}
