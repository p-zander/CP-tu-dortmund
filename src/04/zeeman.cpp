#include <Eigen/Core>
#include <sstream>
#include <utility>
#include <cmath>
#include <functional>
#include <vector>
#include <fstream>

#include "./zeeman.h"
#include "./rungekutta.h"

using Eigen::VectorXd;
using std::endl;
using std::vector;
using std::ofstream;
using std::stringstream;

const double R = 0.05;     // m
const double B = 0.1;      // m
const double A = 0.2;      // m
const double K = 0.25;     // Kg/s²
const double I = 6250e-7;  // Kg*m²
const double G = 2e-4;     // Kg*m²/s

void ab() {
    ofstream file("1_ab.txt");

    file << "# theta \t V(18,1) \t V(16,0) \t V(16,2) \t "
            "F(18,1) ana \t F(18,1) num" << endl;

    for (int i = 0; i < 360; ++i) {
        double theta = M_PI * i / 180;
        file <<        theta         << "\t" << V(theta, 0.18, 0.01) << "\t"
             <<  V(theta, 0.16, 0)   << "\t" << V(theta, 0.16, 0.02) << "\t"
             << F(theta, 0.18, 0.01) << "\t" << -deriveV(theta, 0.18, 0.01)
             << endl;
    }
}

void b(double x, double y, double theta0, double w0 = 0, double tN = 40, size_t steps = 100) {
  vector<double> thetas, ws;

  std::tie(thetas, ws) = RK4_newton(theta0, w0, 0, tN, steps,
                              [x, y](double theta, double w, double) {
                                return 1/I * (-G * w + F(theta, x, y));
                              });

  double h = tN / steps;

  stringstream s;
  s << "1_b_" << x*100 << "_" << y*100 << "_" << theta0 << ".txt";
  ofstream file(s.str());

  for(int i=0; i<steps; i++) {
    file << h*i << "\t" << thetas[i] << "\t" << ws[i] << endl;
  }
}

void c() {
  double x = 0.16;
  double theta0 = 3.88203;
  double w0 = 0;
  size_t steps = 200;

  vector<double> thetas_1, ws_1, thetas_2, ws_2;

  std::tie(thetas_1, ws_1) = RK4_newton(theta0, w0, 0, 100, steps,
                            [x](double theta, double w, double t) {
                                return 1/I * (-G * w + F(theta, x, 0.02*t/100));
                            });
  std::tie(thetas_2, ws_2) = RK4_newton(thetas_1.back(), w0, 0, 100, steps,
                            [x](double theta, double w, double t) {
                                return 1/I * (-G * w + F(theta, x, 0.02*(1-t/100)));
                            });

  ofstream file("1_c.txt");
  double h = 100./steps;

  for (int i=0; i < steps; i++)
    file <<   h*i   << "\t" << thetas_1[i] << "\t" << ws_1[i] << endl;
  for(int i=0; i < steps; i++)
    file << h*i+100 << "\t" << thetas_2[i] << "\t" << ws_2[i] << endl;
}

int main() {
  ab();

  b(0.16, 0, 2);
  b(0.16, 0, 4);
  b(0.16, 0.02, 2);
  b(0.16, 0.02, 4);

  c();

  return 0;
}
