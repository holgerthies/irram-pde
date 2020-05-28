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
  std::function<REAL(const std::array<unsigned int, 2>&, const std::array<REAL,2>&)> fx = [] (auto index, auto x) {
                                                                                            if(index[1] > 0) return REAL(0);
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
  
  std::function<REAL(const std::array<unsigned int, 2>&, const std::array<REAL,2>&)> fy = [] (auto index, auto x) {
                                                                                            if(index[0] > 0) return REAL(0);
                                                                                            switch(index[1] % 4){
                                                                                            case 0:
                                                                                              return sin(x[1]);
                                                                                            case 1:
                                                                                              return cos(x[1]);
                                                                                            case 2:
                                                                                              return -sin(x[1]);
                                                                                            case 3:
                                                                                              return -cos(x[1]);
                                                                                            }
                                                                                            return REAL(0);
                                                                                          };
  std::shared_ptr<Cinfinity<1,REAL>> sine = std::make_shared<Cinfinity<1,REAL>>(f, [] (auto index, auto x, auto eps) {return 1;} );
  std::shared_ptr<Cinfinity<2,REAL>> sinex = std::make_shared<Cinfinity<2,REAL>>(fx, [] (auto index, auto x, auto eps) {if(index[1] > 0) return 0; else return 1;} );
  std::shared_ptr<Cinfinity<2,REAL>> siney = std::make_shared<Cinfinity<2,REAL>>(fy, [] (auto index, auto x, auto eps) {if(index[0] > 0) return 0; else return 1;} );
  auto cosiney = std::make_shared<Cinfinity<2,REAL>>(siney->derive({0,1}));
  auto sxcy = sinex*cosiney;
  // cout << (*sxcy)({0.2,0.3}) << std::endl;
  // cout << std::endl;
  // auto cosine = std::make_shared<Cinfinity<1,REAL>>(sine->derive({1}));
  std::array<REAL,2> center = {{0,0}};
  // std::array<REAL,1> center = {{0}};
  PS_ptr<2,REAL> sp = std::make_shared<Powerseries<2,REAL>>(sinex,center,REAL(1));
  cout <<"test:" << sp({0.1,0.1}) << std::endl;
  PS_ptr<2,REAL> cp = std::make_shared<Powerseries<2,REAL>>(cosiney,center,1);
  auto test = sp->derive({1,0});
  // cout << (sp)({0,0}) << std::endl;
  cout << (cp)({0,0.1}) << std::endl;
  // auto  sp2  = sp*(sp*cp+cp)+sp;
  MVPowerseries<2,2,2,REAL> fun({{{sp,cp},{cp, sp}}});
  MVPowerseries<2,2,1,REAL> v({{{sp},{cp}}});
  print(fun({"0.1", "0.2"}));
  // auto dv = fun.partial_derivative(1);
  // print(dv({"0.1", "0.2"}));
  // cout << std::endl;
  std::array<MVPowerseries<2,2,2,REAL>,2> farr = {{fun,fun}};
  auto sol = Pde_local_solution<2,2,REAL>(farr, v, center);
  print(sol({0.2}));
  
}
