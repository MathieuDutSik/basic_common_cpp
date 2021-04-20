#include "COMB_Stor.h"




int main()
{
  size_t n=10;
  size_t n_iter = 20;
  vector_face vf(n);

  for (size_t iter=0; iter<n_iter; iter++) {
    Face f(n);
    for (size_t i=0; i<n; i++) {
      f[i] = rand() % 2;
    }
    vf.InsertFace(f);
  }

  for (auto & f : vf) {
    std::cerr << " |f|=" << f.size() << " / " << f.count() << "\n";
  }
}

