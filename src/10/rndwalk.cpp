#include <cstdio>
#include <random>
#include <sstream>
#include <vector>


void rnd_walk(double a, size_t n_walks, size_t steps){

  std::mt19937 gen;
  std::uniform_real_distribution<double> dist(-a, std::nextafter(a, a + 1));

  auto rnd = [&dist, &gen]() { return dist(gen); };

  std::vector<std::vector<double>> data(2, std::vector<double>(n_walks, 0));
  
  std::stringstream f1("");
  std::stringstream f2("");
  f1 << "rnd_a_" << a << ".txt";
  f2 << "rnd_a_" << a << "_hist.txt";


  FILE* fout1 = fopen(f1.str().c_str(), "w");
  FILE* fout2 = fopen(f2.str().c_str(), "w");
  fprintf(fout1, "#t\trx1\try1\tr2\n");

  for (size_t i = 1; i <= steps; i++) {
    // #pragma omp parallel for
    for (size_t j = 0; j < n_walks; j++) {
      data[0][j] += rnd();
      data[1][j] += rnd();
    }
    if (a == 1.0 && (i == 10 or i == 100 or i == 1000)){
      for(auto it1 = data[0].begin(),it2=data[1].begin(); it1 != data[0].end();it1++,it2++){
        fprintf(fout2, "%f\t%f\t", *it1, *it2);
      }
      fprintf(fout2, "\n");
    }

    double rx = std::accumulate(data[0].begin(), data[0].end(), 0.0) / n_walks;
    double ry = std::accumulate(data[1].begin(), data[1].end(), 0.0) / n_walks ;
    double r2 = std::inner_product(data[0].begin(), data[0].end(),
                                     data[0].begin(), 0.0) / n_walks +
                  std::inner_product(data[1].begin(), data[1].end(),
                                     data[1].begin(), 0.0) / n_walks;

    fprintf(fout1, "%f\t%f\t%f\n", rx, ry,r2);     

  }
  fclose(fout1);
  fclose(fout2);


}


int main() {
  const double a[4] = {1.0, 1.0e-1, 1.0e-2  , 10.0};
  const size_t n_walks = 100000;
  const size_t steps = 1000;

  for (int i = 0; i < 1; i++){
  rnd_walk(a[i],  n_walks, steps);
  }
  return 0;
}