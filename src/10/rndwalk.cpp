#include <cstdio>
#include <random>
#include <sstream>
#include <vector>


void rnd_walk(double a, size_t n_walks, size_t steps, std::string filename){

  std::mt19937 gen;
  std::uniform_real_distribution<double> dist(-a, std::nextafter(a, a + 1));

  auto rnd = [&dist, &gen]() { return dist(gen); };

  std::vector<std::vector<double>> data(2, std::vector<double>(steps, 0));
  
  FILE* fout = fopen(filename.c_str(), "w");
  fprintf(fout, "#t\trx1\try1\t....\n");

  for (size_t i = 1; i <= n_walks; i++) {
    // #pragma omp parallel for
    for (size_t j = 0; j < steps; j++) {
      data[0][j] += rnd();
      data[1][j] += rnd();
    }
      for(auto it1 = data[0].begin(),it2=data[1].begin(); it1 != data[0].end();it1++,it2++){
      fprintf(fout, "%f\t%f\t", *it1, *it2);
      }
      fprintf(fout, "\n");
  }
  fclose(fout);


}


int main() {
  const double a[4] = {1.0, 1.0e-1, 1.0e-2  , 10.0};
  const size_t n_walks = 1000;
  const size_t steps = 1000;

  for (int i = 0; i < 1; i++){
  std::stringstream filename("");
  filename << "rnd_a_" << a[i] << ".txt";
  rnd_walk(a[i],  n_walks, steps, filename.str());
  }
  return 0;
}