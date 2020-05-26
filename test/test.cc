#include <iRRAM.h>
#include "cinfinity.h"
#include "matrix.h"
#include "powerseries.h"
#include "pde.h"
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
  FactorCache::init();
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
  std::shared_ptr<Cinfinity<1,REAL>> sine = std::make_shared<Cinfinity<1,REAL>>(f, [] (auto index, auto x, auto eps) {return 1;} );
  auto cosine = std::make_shared<Cinfinity<1,REAL>>(sine->derive({101}));
  std::array<REAL,1> center = {{0}};
  PS_ptr<1,REAL> sp = std::make_shared<Powerseries<1,REAL>>(sine,center,REAL(1));
  PS_ptr<1,REAL> cp = std::make_shared<Powerseries<1,REAL>>(cosine,center,1);
  cout << (sp)({0.2}) << std::endl;
  cout << (cp)({0.2}) << std::endl;
  auto  sp2  = (sp*cp)*sp+cp*sp;
  
  cout << (sp2)({0.2}) << std::endl;
  MVPowerseries<1,2,2,REAL> fun({{{sp,cp},{cp, sp}}});
  MVPowerseries<1,2,1,REAL> v({{{sp},{cp}}});
  print(fun({0.2}));
  auto dv = fun.partial_derivative(1);
  print(dv({0.2}));
  cout << std::endl;
  std::array<MVPowerseries<1,2,2,REAL>,1> farr = {{fun}};
  auto sol = Pde_local_solution<2,1,REAL>(farr, v, center);
  print(sol({0.2}));
  
}
