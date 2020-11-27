#include <iRRAM.h>
#include "cinfinity.h"
#include "operator.h"
#include "funops.h"
// #include "matrix.h"
#include "powerseries.h"
#include "mfun.h"
#include "polynomial.h"
#include "diffop.h"
#include "pde.h"
#include "pde_examples.h"
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

// template<unsigned int e, unsigned int d>
// std::vector<REAL> solve(Pde_local_solution<e,d,REAL>& system, const REAL& x){
//   std::vector<REAL> sols(e);
//   auto y = system({x});
//   for(int i=0; i<e;i++){
//     sols[i] = y(i,0);
//   }
//   return sols;
// }

template<unsigned int d>
MVFunction<d,1,1,REAL> get_heat_v(){
  return MVFunction<d,1,1,REAL>();
}
template <>
MVFunction<1,1,1,REAL> get_heat_v<1>(){
    return P2M<1,1,1,REAL>({{{"3x+-1x^2"}}}, {'x'});
}
template <>
MVFunction<2,1,1,REAL> get_heat_v<2>(){
    return P2M<2,1,1,REAL>({{{"x^3+y^2+xy+2"}}}, {'x','y'});
}
template <>
MVFunction<3,1,1,REAL> get_heat_v<3>(){
    return P2M<3,1,1,REAL>({{{"x^3+y^2+xyz^3+2+z^7"}}}, {'x','y','z'});
}

template<unsigned int d>
std::vector<REAL> get_heat_solution(const std::vector<REAL>& x, const REAL& t){
  REAL alpha=REAL(1)/REAL(5);
  auto v = get_heat_v<d>();
  auto diff = heat_equation<d>(alpha);
  vector<REAL, d> xv;
  for(int i=0; i<d;i++) xv[i] = x[i];
  auto sol = solve_pde(diff, v, xv, 10,2);
  std::vector<REAL> ys;
  for(auto& s : sol){
    ys.push_back(s.sum({t}));
  }
  return ys;
}

void compute(){
  FactorCache::init();
  auto D_le = linear_elasticity(3,2,5);
   static int iteration_counter = 0;

   int dimension=0, system,max_iter,prec;
   struct rusage usage;
   struct timeval start, end;
   iRRAM::cout << "choose system" << std::endl;
   iRRAM::cin >>  system;
   while(dimension == 0){
     iRRAM::cout << "choose system dimension" << std::endl;
     iRRAM::cin >> dimension;
     if(dimension != 1 && dimension != 2 && dimension != 3){
       iRRAM::cout << "invalid dimension" << std::endl;
       dimension = 0;
     }
   }
  
   std::vector<REAL> xs;
   iRRAM::cout << "choose x" << std::endl;
   for(int i=0; i<dimension; i++){
     REAL xi;
     iRRAM::cin >> xi;
     xs.push_back(xi);
   }
   REAL t;
   iRRAM::cout << "choose t" << std::endl;
   iRRAM::cin >>  t;

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
     ys = get_heat_solution<1>(xs, t); 
   }
   if(dimension == 2){
     ys = get_heat_solution<2>(xs, t); 
   }
   if(dimension == 3){
     ys = get_heat_solution<3>(xs, t); 
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
  
  
   std::cout << " dimension:" << dimension;
   std::cout << " system:" << system;
   std::cout << " time:" << iteration_time;
   std::cout << " iteration:" << iteration_counter;
   std::cout << " precision:" << ACTUAL_STACK.actual_prec;
   std::cout << " error_mantissa:"  << error.mantissa;
   std::cout << " error_exponent:" << error.exponent;
   std::cout << " normalized:" << error_exp_normalized;
  
   std::cout << std::endl;
  
   if(prec <= 0)
     iRRAM::cout << (REAL(0) == REAL(0)) << std::endl;
  
}
