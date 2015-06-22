#include <algorithm>
#include <cstdio>
#include <numeric>
#include <random>
#include <vector>

int main() {
  const double a = 1.0e-2;
  const size_t N = 1000;
  const size_t n = 100000;
  std::mt19937 gen;
  std::uniform_real_distribution<double> dist(-a, std::nextafter(a, a + 1));

  auto rnd = [&dist, &gen]() { return dist(gen); };

  std::vector<std::vector<double>> data(2, std::vector<double>(N, 0));
  std::vector<double> r(2 * n, 0);
  std::vector<double> r2(n, 0);

  for (size_t i = 2; i < 2 * n; i+=2) {
    for (size_t j = 0; j < N; ++j) {
      data[0][j] += rnd();
      data[1][j] += rnd();
    }
    r[i] = std::accumulate(data[0].begin(), data[0].end(), 0.0) / (N * a);
    r[i + 1] = std::accumulate(data[1].begin(), data[1].end(), 0.0) / (N * a);


    r2[i / 2] = std::inner_product(data[0].begin(), data[0].end(),
                                   data[0].begin(), 0.0) / (N * a * a) +
                std::inner_product(data[1].begin(), data[1].end(),
                                   data[1].begin(), 0.0) / (N * a * a);
  }

  FILE* fout = fopen("rnd.txt", "w");
  fprintf(fout, "#<x>/a\t<y>/a\t<r2>/a\n");
  for (auto it2 = r2.begin(), it = r.begin(); it != r.end(); it2++,it += 2) {
    fprintf(fout, "%f\t%f\t%f\n", *it, *(it + 1), *it2);
  }

  fclose(fout);

  return 0;
}