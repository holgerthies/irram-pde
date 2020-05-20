#include <iRRAM.h>
#include "cinfinity.h"
#include "matrix.h"
using namespace iRRAM;
template<unsigned int m, unsigned int n>
void print(const Matrix<m,n,REAL>& M){
  for(int i=0; i<m; i++){
    for(int j=0;j<n;j++){
      cout << M(i,j) << " ";
    }
    cout << std::endl;
  }
}

void compute(){
  std::function<REAL(const std::array<unsigned int, 1>&, const std::array<REAL,1>&)> f = [] (auto index, auto x) {
                                                                                           switch(index[0] % 4){
                                                                                           case 0:
                                                                                             return sin(x[0]);
                                                                                           case 1:
                                                                                             return cos(x[0]);
                                                                                           case 2:
                                                                                             return -sin(x[0]);
                                                                                           case 3:
                                                                                             return -cos(x[0]);
                                                                                           }
                                                                                           return REAL(0);
                                                                                         };
  Cinfinity<1,REAL> sine(f, [] (auto index, auto x, auto eps) {return 1;} );
  Cinfinity<1,REAL> cosine = sine.derive({101});
  MVFunction<1,2,2,REAL> fun({{{sine,cosine},{cosine, sine}}});
  MVFunction<1,2,1,REAL> v({{{sine},{cosine}}});
  print((v+v)({0.2}));
  cout << std::endl;
  print((fun*v)({0.2}));
  cout << std::endl;
  print((fun({0.2})*v({0.2})));
  auto g = (sine * (cosine + sine)*cosine);
  auto dg = g.derive({50});
  cout << dg({0.2}) << std::endl;
  cout << (cosine)({0.2}) << std::endl;
}
