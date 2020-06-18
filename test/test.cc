#include <iRRAM.h>
#include "mfun.h"
#include "matrix.h"
#include "powerseries.h"
#include "pde.h"
#include "polynomial.h"
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
  std::shared_ptr<Cinfinity<2,REAL>> sinex = std::make_shared<Cinfinity<2,REAL>>(fx, [] (auto index, auto x, auto eps) {  if(index[1] > 0) return 0; else return 1;} );
  std::shared_ptr<Cinfinity<2,REAL>> siney = std::make_shared<Cinfinity<2,REAL>>(fy, [] (auto index, auto x, auto eps) { if(index[0] > 0) return 0; else return 1;} );
  auto cosiney = std::make_shared<Cinfinity<2,REAL>>(siney->derive({0,1}));
  std::array<REAL,2> center = {{0.5,0.75}};
  auto sxcy = toPowerseries<2,REAL>((sinex+cosiney+sinex*cosiney-sinex)*cosiney, center, 1);
  //auto x2 = std::make_shared<Polynomial<1,REAL>>(monomial<REAL>(2,5));
  REAL r=1,c=1;
  auto f0 =(1/r)*monomial<REAL>(0,2);
  auto x = monomial<REAL>(1,2);
  PS_ptr<1,REAL> f12 = std::make_shared<Powerseries<1,REAL>>(f0);
  auto f2 = (r*c*c)*monomial<REAL>(0,2);
  PS_ptr<1,REAL> f21 = std::make_shared<Powerseries<1,REAL>>(f2);
  PS_ptr<1,REAL> zero = std::make_shared<Powerseries<1,REAL>>();
  PS_ptr<1,REAL> lf = std::make_shared<Powerseries<1,REAL>>(x);
  std::array<REAL,1> center1d = {{0}};
  auto F1 = P2M<2,3,3,REAL>({
      {{"0", "0", "1"},
       {"0","0", "0"},
       {"1","0","0"}
      }}, {'x', 'y'});
  auto F2 = P2M<2,3,3,REAL>({
      {{"0", "0", "0"},
       {"0","0", "1"},
       {"0","1","0"}
      }}, {'x', 'y'});
  auto v = P2M<2,3,1,REAL>({{{"2xy+1"}, {"1x^2+0.5"}, {"4y+3"}}}, {'x', 'y'});
  print(F1({0.2,0}));
  print(F2({0.2,0.2}));
  print(v({0.2,0.2}));
  //auto p1 = Polynomial<2,REAL>("7.4+53x^3+23x^2y^1+77y^2x^9+4x^2y^4+6xy+8", {'x', 'y'}, 1);
  //p1.print({'x', 'y'});
//   MVPowerseries<1,3,3,REAL> f_1({{{zero,zero,f12},{zero,zero,zero},{f21,zero, zero}}});
//   MVPowerseries<1,3,3,REAL> f_2({{{zero,zero,zero},{zero,zero,f12},{zero,f21, zero}}});
//   //auto f12 = (1/r)*constf;
// // auto x22 = from1d<2,REAL>(x2, 0);
// // auto y2 = from1d<2,REAL>(x, 1);
// // auto xy = x22*y2;
// // cout << (*sxcy)({0.2,0.3}) << std::endl;
// // cout << std::endl;
//    auto cosine = std::make_shared<Cinfinity<1,REAL>>(sine->derive({1}));
//    PS_ptr<1,REAL> sn = std::make_shared<Powerseries<1,REAL>>(sine,center1d,REAL(1));
// // auto cosine = std::make_shared<Cinfinity<1,REAL>>(sine->derive({1}));
//    PS_ptr<1,REAL> cn = std::make_shared<Powerseries<1,REAL>>(cosine,center1d,REAL(1));
// // PS_ptr<2,REAL> sp = std::make_shared<Powerseries<2,REAL>>(sinex,center,REAL(1));
// // //cout <<"test:" << sp({0.1,0.1}) << std::endl;
// // PS_ptr<2,REAL> cp = std::make_shared<Powerseries<2,REAL>>(cosiney,center,1);
// // // cout << (sp)({0,0}) << std::endl;
// // // auto  sp2  = sp*(sp*cp+cp)+sp;
// // MVPowerseries<1,2,2,REAL> fun1d({{{s,c},{c, s}}});
//    MVPowerseries<1,3,1,REAL> v1d({{{lf},{cn}, {cn}}});
     std::array<MVPowerseries<2,3,3,REAL>,2> farr1d = {{F1,F2}};
     auto sol = Pde_local_solution<3,2,REAL>(farr1d, v, {{0,0}});
     for(int i=0; i<1000;i++){
       print(sol({REAL(i)/REAL(100)}));
       cout << std::endl;
     }

//print(fun1d({0.2}));
// MVPowerseries<2,2,2,REAL> fun({{{sp,cp},{cp, sp}}});
// MVPowerseries<2,2,1,REAL> v({{{sp},{cp}}});
// MVPowerseries<2,2,2,REAL> fun({{{sp,cp},{cp, sp}}});
// MVPowerseries<2,2,1,REAL> v({{{sp},{cp}}});
// // print(fun({"0", "0"}));
// // auto dv = fun.partial_derivative(1);
// // print(dv({"0.1", "0.2"}));
// // cout << std::endl;
// std::array<MVPowerseries<2,2,2,REAL>,2> farr = {{fun,fun}};
// auto sol = Pde_local_solution<2,2,REAL>(farr, v, center);
// print(sol({1}));
  
}
