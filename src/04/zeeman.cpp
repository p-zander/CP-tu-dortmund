#include <iostream>
#include <cmath>
#include <functional>

using std::cout;
using std::endl;

const double R = 0.05;     // m
const double B = 0.1;      // m
const double A = 0.2;      // m
const double K = 0.25;     // Kg/s²
// const double I = 6250e-7;  // Kg*m²
// const double G = 2e-4;     // Kg*m²/s

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


int main() { 

	const double x = 0.18;
	const double y = 0.01;

	for (int i = 0; i < 360 ; ++i)
	{
          double theta = M_PI * i / 180;
          cout << i << "\t" << V(theta, x, y) << "\t" << F(theta, x, y) << "\t" << -deriveV(theta, x, y) << endl;
        }

	return 0; 
}