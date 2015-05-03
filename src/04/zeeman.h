#include <functional>
#include <cmath>

extern const double R;
extern const double B;
extern const double A;
extern const double K;
extern const double I;
extern const double G;

double V(double theta, double x, double y) {
  double d1 = sqrt(A * A + 2 * A * R * cos(theta) + R * R);
  double d2 = sqrt(R * R + x * x + y * y - 2 * R * (x * cos(theta) + y * sin(theta)));
  return .5 * K * pow(d1 - B, 2) + .5 * K * pow(d2 - B, 2);
}

double F(double theta, double x, double y) {
  double d1 = sqrt(A * A + 2 * A * R * cos(theta) + R * R);
  double d2 = sqrt(R * R + x * x + y * y - 2 * R * (x * cos(theta) + y * sin(theta)));
  double z1 = -A * R * sin(theta);  // numerator of derivative of d1
  double z2 = -R * (y * cos(theta) - x * sin(theta));  // numerator of derivative of d2
  return -K * (z1 - B * z1 / d1) - K * (z2 - B * z2 / d2);
}

double deriveV(double theta, double x, double y, double h = 1e-3) {
  auto V2 = std::bind(V, std::placeholders::_1, x, y);
  return (V2(theta - 2 * h) - 8 * V2(theta - h) + 8 * V2(theta + h) - V2(theta + 2 * h)) / (12 * h);
}
