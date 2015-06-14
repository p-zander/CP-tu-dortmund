#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/src/KroneckerProduct/KroneckerTensorProduct.h>
#include <iostream>
#include <cmath>
#include <utility>

using Eigen::Matrix;

template <typename T, size_t N>
std::pair<T, Matrix<T, N, 1>> powerEig(const Matrix<T, N, N> &A, Matrix<T, N, 1> guess, double eps = 1e-6) {
    Matrix<T, N, 1> copy;
    unsigned int steps = 0;

    guess.normalize();

    do {
        steps++;
        copy = guess;
        guess = (A * guess).normalized();
    } while ((copy.cwiseAbs() - guess.cwiseAbs()).unaryExpr([eps](T x) { return x * x > eps * eps; }).any());

    std::cout << "steps needed: " << steps << std::endl;
    return std::make_pair(guess.dot(A * guess), std::move(guess));
}

template <typename T, size_t N>
std::pair<std::vector<T>, std::vector<Matrix<T, N, 1>>> powerEigs(const Matrix<T, N, N> &A, Matrix<T, N, 1> guess,
                                                                  size_t k = 0, double eps = 1e-6) {
    if (k == 0) k = N;

    std::vector<T> vals(N);
    std::vector<Matrix<T, N, 1>> vecs(N);

    Matrix<T, N, N> B(A);
    Matrix<T, Eigen::Dynamic, Eigen::Dynamic> C;

    for (int i = 0; i < k; ++i) {
        std::tie(vals[i], vecs[i]) = powerEig<T, N>(B, guess, eps);
        C = kroneckerProduct(vecs[i], vecs[i]).eval();
        C.resize(N, N);
        B = B - vals[i] * C;
    }

    return std::make_pair(std::move(vals), std::move(vecs));
}

template <typename T> inline T off(const Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &A) {
    Matrix<T, Eigen::Dynamic, Eigen::Dynamic> diag(A.diagonal().asDiagonal());
    return (A - diag).cwiseAbs2().sum();
}

template <typename T, size_t N> void jacobiRotation(Matrix<T, N, N> &A, double eps = 1e-16) {
    Matrix<T, N, N> diag(A.diagonal().asDiagonal());
    Matrix<T, N, N> P;

    double theta, t, c, s;
    size_t p, q;
    unsigned int steps = 0;

    while (off<double>(A) > eps && steps < 500) {
        steps++;
        diag = A.diagonal().asDiagonal();
        (A - diag).cwiseAbs().maxCoeff(&p, &q);
        if (p > q) std::swap(p, q);
        theta = (A(q, q) - A(p, p)) * 1. / (2 * A(p, q));
        t = ((theta > 0) - (theta < 0)) * 1. / (abs(theta) + sqrt(theta * theta + 1));
        c = 1. / sqrt(1 + t * t);
        s = t * c;
        P = Eigen::Matrix<T, N, N>::Identity(N, N);
        P(p, p) = P(q, q) = c;
        P(p, q) = s;
        P(q, p) = -s;
        A = P.transpose() * A * P;
    }

    if (steps == 500)
        std::cerr << std::endl << "Number of steps reached upper limit. Estimate may not be accurate!" << std::endl;
}

int main() {
    Eigen::Matrix4d A;
    A << 1, 2, 2, 4, 2, 5, 2, 1, 2, 2, 6, 3, 4, 1, 3, 4;

    std::cout << "A = " << std::endl << A << std::endl;

    std::cout << std::endl << "A.eigenvalues() = " << std::endl << A.eigenvalues() << std::endl;

    Eigen::Vector4d guess = Eigen::Vector4d::LinSpaced(A.cols(), -1, 1);

    std::cout << std::endl;
    std::cout << "First eigenvalue from powerEig: " << powerEig<double, 4>(A, guess).first << std::endl << std::endl;

    std::vector<double> vals;
    std::tie(vals, std::ignore) = powerEigs<double, 4>(A, guess);

    std::cout << std::endl << "All eigenvalues from powerEigs: " << std::endl;
    for (const double &x : vals) std::cout << x << std::endl;

    std::cout << std::endl << "off(A) = " << off<double>(A) << std::endl;

    jacobiRotation<double, 4>(A);

    std::cout << std::endl << "Jacobi rotated matrix A:" << std::endl << A << std::endl;
    std::cout << std::endl << "Eigenvalues are: " << std::endl << A.diagonal() << std::endl;

    return 0;
}
