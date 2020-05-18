#include <iRRAM.h>
#include "cinfinity.h"
using namespace iRRAM;
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
  auto g = (sine * (cosine + sine)*cosine);
  auto dg = g.derive({50});
  cout << dg({0.2}) << std::endl;
  cout << (cosine)({0.2}) << std::endl;
}
