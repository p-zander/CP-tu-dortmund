#include <Eigen/Core>
#include <sstream>
#include <utility>
#include <cmath>
#include <functional>
#include <vector>
#include <fstream>
#include <array>
#include <string>

#include "zeeman.h"
#include "rungekutta.h"

using Eigen::VectorXd;
using std::endl;
using std::vector;
using std::ofstream;
using std::stringstream;
using std::array;

const double R = 0.05;     // m
const double B = 0.1;      // m
const double A = 0.2;      // m
const double K = 0.25;     // Kg/s²
const double I = 6250e-7;  // Kg*m²
const double G = 2e-4;     // Kg*m²/s

void a() {
  array<double, 6> x{{0.218, 0.219, 0.21925, 0.2195, 0.22, 0.2225}};
  double theta0 = 0;
  double w0 = 0;
  size_t steps = 200;
  double tN = 100;

  vector<vector<double>> thetas, ws;

  for (const double& xi : x) {
    vector<double> thetas_i, ws_i;
    std::tie(thetas_i, ws_i) = RK4_newton(
        theta0, w0, 0, tN, steps, [xi](double theta, double w, double t) {
          return 1 / I *
                 (-G * w + F(theta, xi, 0.025 * sin(2 * M_PI * t / 10)));
        });
    thetas.push_back(thetas_i);
    ws.push_back(ws_i);
  }
  ofstream file("2a.txt");
  file << "t \t theta_i \t omega_i, \t i= 1...6 " << endl;
  for (size_t i = 0; i < steps; ++i) {
    file << i* tN / steps << "\t";
    for (size_t j = 0; j < x.size(); j++) {
      file << thetas[j][i] << "\t" << ws[j][i] << "\t";
    }
    file << endl;
  }
  file.close();
}

void bc(const double x0, const double xN, const size_t xsteps,
        const std::string filename) {
  double theta0 = 0;
  double w0 = 0;
  double tN = 1000;
  double xi = x0;
  size_t steps = 2000;
  double h = (xN - x0) / xsteps;
  vector<vector<double>> thetas, ws;

  for (size_t i = 0; i < xsteps; i++) {
    vector<double> thetas_i, ws_i;
    std::tie(thetas_i, ws_i) = RK4_newton(
        theta0, w0, 0, tN, steps, [xi](double theta, double w, double t) {
          return 1 / I *
                 (-G * w + F(theta, xi, 0.025 * sin(2 * M_PI * t / 10)));
        });
    thetas.push_back(thetas_i);
    ws.push_back(ws_i);
    xi += h;
  }
  ofstream file(filename);
  file << "x \t theta \t omega" << endl;
  for (size_t j = 0; j < xsteps; j++) {
    for (size_t i = 100.*steps/tN ; i < steps; i += 10.*steps/tN) {
      file << x0 + j* h << "\t" << thetas[j][i] << "\t" << ws[j][i] << endl;
    }
  }
  file.close();
}

int main() {
  a();
  bc(0.22, 0.22, 1, "2b.txt");
  bc(0.218, 0.22, 200, "2c.txt");

  return 0;
}