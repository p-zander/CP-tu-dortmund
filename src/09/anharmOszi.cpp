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
    H += (0.2 * VectorXd::LinSpaced(N, -L, L).cwiseAbs2().cwiseAbs2()).asDiagonal();

    Eigen::VectorXd vals(H.eigenvalues().real());
    std::sort(vals.data(), vals.data() + vals.size());  // eigenvalues are unsorted

    std::cout << vals.head(10) << std::endl << std::endl;
}

void c(const size_t N) {
    VectorXi M(VectorXi::LinSpaced(N, -(N - 1) / 2, (N - 1) / 2));

    MatrixXd H = M.cast<double>().unaryExpr([](double m) { return 6 * m * m + 6 * m + 3; }).asDiagonal();

    M = M.segment(1, N - 2);
    H.bottomLeftCorner(N - 2, N - 2) +=
        M.cast<double>().unaryExpr([](double m) { return (4 * m - 2) * sqrt(m * (m - 1)); }).asDiagonal();
    H.topRightCorner(N - 2, N - 2) +=
        M.cast<double>().unaryExpr([](double m) { return (4 * m + 6) * sqrt((m + 1) * (m + 2)); }).asDiagonal();

    M = M.segment(1, N - 4);
    H.bottomLeftCorner(N - 4, N - 4) +=
        M.cast<double>().unaryExpr([](double m) { return sqrt(m * (m - 1) * (m - 2) * (m - 3)); }).asDiagonal();
    H.topRightCorner(N - 4, N - 4) +=
        M.cast<double>().unaryExpr([](double m) { return sqrt((m + 1) * (m + 2) * (m + 3) * (m + 4)); }).asDiagonal();

    H *= 1. / 4;

    Eigen::VectorXd vals(H.eigenvalues().real());
    std::sort(vals.data(), vals.data() + vals.size());  // eigenvalues are unsorted

    std::cout << vals.head(10) << std::endl;
}

int main() {
    const double L{10};  // computation window is [-L, +L] in both directions
    const double d_l{0.1};

    const int N(1 + static_cast<int>(2 * L / d_l + 0.5));

    b(N, L);
    c(50);

    return 0;
}
