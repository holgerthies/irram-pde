using REAL = double;
using INTEGER = long;
#include <iostream>
#include <limits>

using namespace std;
namespace iRRAM{
double power(const double& a, const int& b) {return pow(a,b);}
}
double maximum(const double& a, const double& b) { return max(a,b); }
#include "combinatorics_double.h"
#include "cinfinity.h"
#include "operator.h"
#include "funops.h"
// #include "matrix.h"
#include "powerseries_double.h"
#include "mfun.h"
#include "polynomial.h"
#include "diffop.h"
#include "pde.h"
#include "pde_examples.h"
using namespace iRRAM;

template<unsigned int m, unsigned int n>
void print(const Matrix<m,n,double>& M){
  for(int i=0; i<m; i++){
    for(int j=0;j<n;j++){
      std::cout << M(i,j) << " ";
    }
    std::cout << std::endl;
  }
}
typedef std::numeric_limits< double > dbl;
int main(){
  FactorCache::init();  
  // REALMATRIX A = zeroes(2,2), B=zeroes(2,2);
  // A(0,0) = 6;
  // A(1,1) = 1/A(0,0);
  // B(0,1) = 2;
  // B(1,0) = B(0,1);

  // cout << to_string_double(symm_eig(A)) << std::endl;
  // cout << to_string_double(symm_eig(B)) << std::endl;
  // auto s = pencil_decomposition(A, B, 2, 2);
  // cout << to_string_double(s.A.eigenvalues) << std::endl;
  // cout << to_string_double(s.B.eigenvalues) << std::endl;
  // cout << to_string_double(s.B.eigenvectors) << std::endl;



  auto D_le = linear_elasticity(3,2,5);
   int dimension=0, system;
   struct rusage usage;
   struct timeval start, end;
   auto f = sine<3>(1,0.5);
   auto f2 = sine<2>(1,0.5);
   auto g = sine<3>(2,0.3);
   auto g2 = sine<2>(0,1);
   auto h = sine<3>(0,0.5);
   auto v = R2M<3,9,1,REAL>({0.1,0.,0.3,0.4,0.5,0.6,0.7,0.8,0.9});
   v(0,0) = f;
   v(3,0) = f*g;
   v(4,0) = h+g+f;
   auto v2 = R2M<2,3,1,REAL>({REAL(0.1),REAL(0.2),REAL(0.3)});
   v2(0,0) = f2;
   v2(1,0) = f2+g2;
   auto elasticity = linear_elasticity(1,2,3);
   auto elasticity_const = linear_elasticity_const(1,2,3);
   auto acoustics = acoustics_2d<REAL>(1,2);
   auto acoustics_const = acoustics_2d_const(1,2);
  // auto v2 = v;
  // auto Di = system1c;
  // for(int i=0; i<100; i++){
  //   cout << "i:" << i << std::endl;
  //   //v2 = system1c(v2);
  //   //print(v2.get_derivative({0,0,0}));
  //   cout << std::endl;
  //   cout <<"Equals:"<< std::endl;
  //   cout << std::endl;
  //   print(Di(v).get_derivative({0,0,0}));
  //   Di = system1c*Di;
  // }
  cout.precision(dbl::max_digits10);
   auto system2 = linear_elasticity_trig();
   cout << "choose example" << std::endl;
   cout << "1) 2D acoustics with constant coefficients" << std::endl;
   cout << "2) 2D acoustics with constant function coefficients" << std::endl;
   cout << "3) 3D linear elasticitiy with constant coefficients" << std::endl;
   cout << "4) 3D linear elasticitiy with constant coefficients (method 2)" << std::endl;
   cout << "5) 3D linear elasticitiy with constant function coefficients" << std::endl;
   cout << "6) 3D linear elasticitiy with trigonometric coefficients" << std::endl;
   cin >>  system;
  
   iRRAM::vector<REAL, 3> x;
   iRRAM::vector<REAL, 2> x2;
   REAL r_max;
   switch(system){
     case 1: 
      dimension = 2;
      r_max = get_sol_r_const(2, 3, 1, 1);
     break;
     case 2: 
      dimension = 2;
      r_max = get_sol_r(2, 3, 1, 1);
     break;
     case 3: 
      dimension = 3;
      r_max = get_sol_r_const(3, 9, 1, 1);
     break;
     case 4: 
      dimension = 3;
      r_max = get_sol_r_const(3, 9, 1, 1);
     break;
     case 5: 
      dimension = 3;
      r_max = get_sol_r(3, 9, 1, 1);
     break;
     case 6: 
      dimension = 3;
      r_max = get_sol_r(3, 9, 1, 1);
     break;
   }
   cout << "choose x "  << std::endl;
   for(int i=0; i<dimension; i++){
     cin >> x[i];
   }
   for(int i=0; i<2; i++){
      x2[i] = x[i];
   }
   REAL t;
   cout << "choose t < " << r_max << std::endl;
   cin >>  t;

   getrusage(RUSAGE_SELF, &usage);
   start = usage.ru_utime;
  
   getrusage(RUSAGE_SELF, &usage);
   start = usage.ru_utime;
   std::vector<REAL> ys;
   
   if(system == 1){
     auto sol = solve_pde_constant_method2(acoustics_const, v2, x2, 1, 1);
     for(auto& s : sol){
       ys.push_back(s.sum({t}));
     }
   }
   if(system == 2){
     auto sol = solve_pde(acoustics, v2, x2, 1, 1);
     for(auto& s : sol){
       ys.push_back(s.sum({t}));
     }
   }
   if(system == 3){
     auto sol = solve_pde_constant(elasticity_const, v, x, 1, 1);
     for(auto& s : sol){
       ys.push_back(s.sum({t}));
     }
   }
   if(system == 4){
     auto sol = solve_pde_constant_method2(elasticity_const, v, x, 1, 1);
     for(auto& s : sol){
       ys.push_back(s.sum({t}));
     }
   }
   if(system == 5){
     auto sol = solve_pde(elasticity, v, x, 1, 1);
     for(auto& s : sol){
       ys.push_back(s.sum({t}));
     }
   }
   if(system == 6){
     auto sol = solve_pde(system2, v, x, 1, 1);
     for(auto& s : sol){
       ys.push_back(s.sum({t}));
     }
   }
   cout << std::endl;
   getrusage(RUSAGE_SELF, &usage);
   end = usage.ru_utime;
   auto iteration_time =  end.tv_sec-start.tv_sec+double(end.tv_usec-start.tv_usec)/1000000;
   for(auto& c : ys){
     std::cout << c << std::endl; 
   }
   std::cout << " dimension:" << dimension;
   std::cout << " system:" << system;
   std::cout << " time:" << iteration_time;
   std::cout << std::endl;
  
}
