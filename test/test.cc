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

template<unsigned int e, unsigned int d>
std::vector<REAL> solve(Pde_local_solution<e,d,REAL>& system, const REAL& x){
  std::vector<REAL> sols(e);
  auto y = system({x});
  for(int i=0; i<e;i++){
    sols[i] = y(i,0);
  }
  return sols;
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
  static int iteration_counter = 0;

  int dimension=0, system,max_iter,prec;
  struct rusage usage;
  struct timeval start, end;
  auto transport_equation = P2MF<1,1,1,REAL>({{{"-3"}}}, {'x'});
  auto l = P2MF<1,1,1,REAL>({{{"x^2+x"}}}, {'x'});
  auto solution = Pde_local_solution<1,1,REAL>({{transport_equation}}, l, {0.3}); 
  std::array<REAL,1> center = {0.3};
  for(int i=0; i<10; i++){
    print(solution({REAL(i)/REAL(30)}));
  }
  auto F1_1df = P2MF<1,2,2,REAL>({
      {{"0", "1"},
       {"1", "0"},
      }}, {'x'});
  auto F1_1d = toPowerseries(F1_1df, {{0}}, 1);

  auto v_1df = P2MF<1,2,1,REAL>({{{"2x+0.02"}, {"1x^2+0.2"}}}, {'x'});
  auto v_1d = toPowerseries(v_1df, {{0}}, 1);

  auto F1_2df = P2MF<2,3,3,REAL>({
      {{"0", "0", "1"},
       {"0","0", "0"},
       {"1","0","0"}
      }}, {'x', 'y'});
  auto F1_2d = toPowerseries(F1_2df, {{0,0}}, 1);

  auto F2_2df = P2MF<2,3,3,REAL>({
      {{"0", "0", "0"},
       {"0","0", "1"},
       {"0","1","0"}
      }}, {'x', 'y'});
  auto F2_2d = toPowerseries(F2_2df, {{0,0}}, 1);

  auto v_2df = P2MF<2,3,1,REAL>({{{"2xy+0.01"}, {"1x^2+0.05"}, {"4y+0.03"}}}, {'x', 'y'});
  auto v_2d = toPowerseries(v_2df, {{0,0}}, 1);

  std::array<MVPowerseries<1,2,2,REAL>,1> farr_1d = {{F1_1d}};
  std::array<MVPowerseries<2,3,3,REAL>,2> farr_2d = {{F1_2d, F2_2d}};

  std::array<MVFunction<1,2,2,REAL>,1> farr_1df = {{F1_1df}};
  std::array<MVFunction<2,3,3,REAL>,2> farr_2df = {{F1_2df, F2_2df}};
  while(dimension == 0){
    iRRAM::cout << "choose system dimension" << std::endl;
    iRRAM::cin >> dimension;
    if(dimension != 1 && dimension != 2){
      iRRAM::cout << "invalid dimension" << std::endl;
      dimension = 0;
    }
  }
  iRRAM::cout << "choose system" << std::endl;
  iRRAM::cin >>  system;
  
  REAL x;
  iRRAM::cout << "choose x" << std::endl;
  iRRAM::cin >>  x;

  iRRAM::cout << "choose precision (or 0 for iteration number)" << std::endl;
  iRRAM::cin >>  prec;
  bool out = true;
  if(prec > 0){
    iRRAM::cout << setRwidth(prec) << std::endl;
  }
  else{
    iRRAM::cout << "choose number of iterations" << std::endl;
    iRRAM::cin >>  max_iter;
    out=false;
  }

  iteration_counter++; 
  if(prec <= 0 && iteration_counter == max_iter) return;
  
  
  getrusage(RUSAGE_SELF, &usage);
  start = usage.ru_utime;
  
  getrusage(RUSAGE_SELF, &usage);
  start = usage.ru_utime;
  std::vector<REAL> ys;
  if(dimension == 1){
    auto sol = Pde_local_solution<2,1,REAL>(farr_1df, v_1df, {{0}});
    ys = solve(sol, x);
  }
  if(dimension == 2){
    auto sol = Pde_local_solution<3,2,REAL>(farr_2d, v_2d, {{0,0}});
    ys = solve(sol, x);
  }
    for(int i = 0; i<ys.size();i++){
      cout << ys[i] << " ";
    }
    cout << std::endl;
  //print(sol({x}));
  getrusage(RUSAGE_SELF, &usage);
  end = usage.ru_utime;
  auto iteration_time =  end.tv_sec-start.tv_sec+double(end.tv_usec-start.tv_usec)/1000000;

  sizetype error;
  sizetype_exact(error);
  for(auto y : ys){
    sizetype curr_error;
    y.geterror(curr_error);
    sizetype_max(error, error, curr_error);
  }
  
  int error_exp_normalized;
  unsigned long mantissa = error.mantissa;
  error_exp_normalized = error.exponent;
  while(mantissa > 1){
    mantissa /= 2;
    error_exp_normalized++;
  }
  //check(sol, A3_sol(5), 50);
  // return;
  
  
  // std::cout << " dimension:" << dimension;
  // std::cout << " system:" << system;
  // std::cout << " time:" << iteration_time;
  // std::cout << " iteration:" << iteration_counter;
  // std::cout << " precision:" << ACTUAL_STACK.actual_prec;
  // std::cout << " error_mantissa:"  << error.mantissa;
  // std::cout << " error_exponent:" << error.exponent;
  // std::cout << " normalized:" << error_exp_normalized;
  std::cout << dimension << "," << -1*error_exp_normalized << "," << iteration_time;
  
  std::cout << std::endl;
  
  if(prec <= 0)
    iRRAM::cout << (REAL(0) == REAL(0)) << std::endl;
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
