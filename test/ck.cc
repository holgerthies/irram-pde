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
#include "matrix_pencil.h"
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
  static int iteration_counter = -1;

   int dimension=0, system,max_iter,prec;
   struct rusage usage;
   struct timeval start, end;
   auto f = sine<3>(1,REAL("0.5"));
   auto f2 = sine<2>(1,REAL("0.5"));
   auto g = sine<3>(2,REAL("0.3"));
   auto g2 = sine<2>(0,1);
   auto h = sine<3>(0,REAL("0.5"));
   auto v = R2M<3,9,1,REAL>({REAL("0.1"),REAL("0.2"),REAL("0.3"),REAL("0.4"),REAL("0.5"),REAL("0.6"),REAL("0.7"),REAL("0.8"),REAL("0.9")});
   v(0,0) = f;
   v(3,0) = f*g;
   v(4,0) = h+g+f;
   auto v2 = R2M<2,3,1,REAL>({REAL("0.1"),REAL("0.2"),REAL("0.3")});
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


   auto system2 = linear_elasticity_trig();
   iRRAM::cout << "choose example" << std::endl;
   iRRAM::cout << "1) 2D acoustics with constant coefficients" << std::endl;
   iRRAM::cout << "2) 2D acoustics with constant function coefficients" << std::endl;
   iRRAM::cout << "3) 3D linear elasticitiy with constant coefficients" << std::endl;
   iRRAM::cout << "4) 3D linear elasticitiy with constant coefficients (method 2)" << std::endl;
   iRRAM::cout << "5) 3D linear elasticitiy with constant function coefficients" << std::endl;
   iRRAM::cout << "6) 3D linear elasticitiy with trigonometric coefficients" << std::endl;
   iRRAM::cin >>  system;
  
   vector<REAL, 3> x;
   vector<REAL, 2> x2;
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
   iRRAM::cout << "choose x "  << std::endl;
   for(int i=0; i<dimension; i++){
     iRRAM::cin >> x[i];
   }
   for(int i=0; i<2; i++){
      x2[i] = x[i];
   }
   REAL t;
   iRRAM::cout << "choose t < " << r_max << std::endl;
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

   sizetype error;
   sizetype_exact(error);
   for(auto y : ys){
     sizetype curr_error;
     y.geterror(curr_error);
     sizetype_max(error, error, curr_error);
     //iRRAM::cout << y << std::endl;
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
   std::cout << " precision:" << actual_stack().actual_prec;
   std::cout << " error_mantissa:"  << error.mantissa;
   std::cout << " error_exponent:" << error.exponent;
   std::cout << " normalized:" << error_exp_normalized;
  
   std::cout << std::endl;
  if(out){
    for(auto c : ys){
      iRRAM::cout << c << std::endl;
    }
  }
   if(prec <= 0)
     if(REAL(0) == REAL(0)) iRRAM::cout << "this will never happen" << std::endl;
  
}
