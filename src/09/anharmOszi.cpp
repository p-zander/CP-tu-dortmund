#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;

void b(const size_t N, const double L) {
  const double d_l = 2 * L / (N - 1);

  MatrixXd H = -2 * MatrixXd::Identity(N, N);

  H.topRightCorner(N - 1, N - 1) += MatrixXd::Identity(N - 1, N - 1);
  H.bottomLeftCorner(N - 1, N - 1) += MatrixXd::Identity(N - 1, N - 1);
  H *= -1. / pow(d_l, 2);

  H += VectorXd::LinSpaced(N, -L, L).cwiseAbs2().asDiagonal();
  H += (0.2 * VectorXd::LinSpaced(N, -L, L).cwiseAbs2().cwiseAbs2())
           .asDiagonal();

  Eigen::VectorXd vals(H.eigenvalues().real());
  std::sort(vals.data(),
            vals.data() + vals.size());  // eigenvalues are unsorted

  std::cout << vals.head(10) << std::endl << std::endl;
}

void c(const size_t N) {
  VectorXi M(VectorXi::LinSpaced(N, 0, N));

  MatrixXd H = M.cast<double>()
                   .unaryExpr([](double m) { return 6 * m * m + 6 * m + 3; })
                   .asDiagonal();

  H.topRightCorner(N - 2, N - 2) +=
      M.segment(2, N - 2)
          .cast<double>()
          .unaryExpr([](double m) { return (4 * m - 2) * sqrt(m * (m - 1)); })
          .asDiagonal();
  H.bottomLeftCorner(N - 2, N - 2) +=
      M.segment(0, N - 2)
          .cast<double>()
          .unaryExpr(
               [](double m) { return (4 * m + 6) * sqrt((m + 1) * (m + 2)); })
          .asDiagonal();

  H.topRightCorner(N - 4, N - 4) +=
      M.segment(4, N - 4)
          .cast<double>()
          .unaryExpr(
               [](double m) { return sqrt(m * (m - 1) * (m - 2) * (m - 3)); })
          .asDiagonal();
  H.bottomLeftCorner(N - 4, N - 4) +=
      M.segment(0, N - 4)
          .cast<double>()
          .unaryExpr([](double m) {
            return sqrt((m + 1) * (m + 2) * (m + 3) * (m + 4));
          })
          .asDiagonal();

  H *= 0.25 * 0.2;

  H += M.cast<double>()
           .unaryExpr([](double m) { return 2 * m + 1; })
           .asDiagonal();
  Eigen::VectorXd vals(H.eigenvalues().real());
  std::sort(vals.data(),
            vals.data() + vals.size());  // eigenvalues are unsorted

  std::cout << vals.head(10) << std::endl << std::endl;
}

int main() {
  const double L{10};  // computation window is [-L, +L] in both directions

  const int N1(1 + static_cast<int>(2 * L / 0.5 + 0.5));
  const int N2(1 + static_cast<int>(2 * L / 0.1 + 0.5));
  const int N3(1 + static_cast<int>(2 * L / 0.01 + 0.5));

  std::cout << "Aufgabenteil b)\n"
            << "L=10, dl = 0.5" << std::endl;
  b(N1, L);
  std::cout << "L=10, dl = 0.1" << std::endl;
  b(N2, L);
  std::cout << "L=10, dl = 0.01" << std::endl;
  b(N3, L);
  std::cout << "Aufgabenteil c)\n"
            << "N = 20 " << std::endl;
  c(20);
  std::cout << "N = 50 " << std::endl;
  c(50);
  std::cout << "N = 250 " << std::endl;
  c(250);
  std::cout
      << "\n \n Es scheint so, als wenn die Methode aus Aufgabenteil c genauer ist, da diese schneller gegen den finalen Wert konvergiert. \n";
  return 0;
}
