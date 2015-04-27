#include <Eigen/Core>
#include <functional>
#include <utility>
#include <vector>
#include <cmath>

using Eigen::VectorXd; using std::pair; using std::vector;

pair<VectorXd, VectorXd>
RK4_newton(VectorXd rn, VectorXd vn, double tn, double h,
            const std::function<VectorXd(VectorXd, double)>& f) {
    assert(rn.rows() == vn.rows());

    VectorXd k0 = h *   vn;
    VectorXd l0 = h * f(rn,         tn);
    VectorXd k1 = h *  (vn + .5*l0);
    VectorXd l1 = h * f(rn + .5*k0, tn + .5*h);
    VectorXd k2 = h *  (vn + .5*l1);
    VectorXd l2 = h * f(rn + .5*k1, tn + .5*h);
    VectorXd k3 = h *  (vn +    l2);
    VectorXd l3 = h * f(rn +    k2, tn +    h);

    return std::make_pair(rn + 1./6*(k0 + 2*k1 + 2*k2 + k3), vn + 1./6*(l0 + 2*l1 + 2*l2 + l3));
}

pair<vector<VectorXd>, vector<VectorXd>>
RK4_newton(VectorXd r0, VectorXd v0, double t0, double tN, int N,
            const std::function<VectorXd(VectorXd, double)>& f) {
    pair<VectorXd, VectorXd> val = std::make_pair(r0, v0);
    vector<VectorXd> res_r {r0};
    vector<VectorXd> res_v {v0};

    double h = (tN-t0)/N;

    for (int i=1; i <=N; i++) {
        val = RK4_newton(std::get<0>(val), std::get<1>(val), h*i, h, f);
        res_r.push_back(std::get<0>(val));
        res_v.push_back(std::get<1>(val));
    }

    return std::make_pair(res_r, res_v);
}

