#ifndef PDE_EXAMPLES_H
#define PDE_EXAMPLES_H
#include "cinfinity.h"
#include "operator.h"
#include "funops.h"
#include "powerseries.h"
#include "mfun.h"
#include "polynomial.h"
#include "diffop.h"
#include "pde.h"
namespace iRRAM{
  // returns sin(a*x_i)
  template<unsigned int d>
  std::shared_ptr<Cinfinity<d,REAL>> sine(const int i, const REAL& a){
    return std::make_shared<Cinfinity<d,REAL>>([i, a] (const Multiindex<d>& m, const vector<REAL,d>& x) {
                                                 for(int k=0; k<d; k++){
                                                   if(k != i && m[k] != 0) return REAL(0);
                                                 }
                                                 auto k=m[i];
                                                 auto fact = power(a, k); 
                                                 switch(k % 4){
                                                 case 0:
                                                   return fact*sin(a*x[i]);
                                                 case 1:
                                                   return fact*cos(a*x[i]);
                                                 case 2:
                                                   return -fact*sin(a*x[i]);
                                                 case 3:default:
                                                   return -fact*cos(a*x[i]);
                                                 }
                                               });
  }
    // 2D acoustics
  template<class T>
  DifferentialOperator<2,3,T> acoustics_2d(const T& r, const T& c){
    std::array<std::array<REAL, 3>,3> M1d({{
            {0, 0, REAL(1)/r},
            {0, 0, 0},
            {r*c*c, 0, 0},
          }});
    MVFunction<2,3,3,T> M1(R2M<2,3,3,T>(M1d));

    std::array<std::array<REAL, 3>,3> M2d({{
            {0, 0, 0},
            {0, 0, T(1)/r},
            {0, r*c*c, 0},
          }});
    MVFunction<2,3,3,T> M2(R2M<2,3,3,T>(M2d));
   
   auto d1 = DifferentialOperator<2,3,T>({1,0},M1);
   auto d2 = DifferentialOperator<2,3,T>({0,1},M2);
   return d1+d2;
  }

    // 2D acoustics
  ConstantDifferentialOperator<2,3,REAL> acoustics_2d_const(const REAL& r, const REAL& c){
    std::array<std::array<REAL, 3>,3> M1d({{
            {0, 0, REAL(1)/r},
            {0, 0, 0},
            {r*c*c, 0, 0},
          }});
    auto M1 = Matrix<3,3,REAL>(M1d);

    auto M2 = Matrix<3,3,REAL>({{
            {0, 0, 0},
            {0, 0, REAL(1)/r},
            {0, r*c*c, 0},
          }});
   
   auto d1 = ConstantDifferentialOperator<2,3,REAL>({1,0},M1);
   auto d2 = ConstantDifferentialOperator<2,3,REAL>({0,1},M2);
   return d1+d2;
  }

  // linear elasticity equation
  DifferentialOperator<3,9,REAL> linear_elasticity(const REAL& r, const REAL& L, const REAL& M){
    std::array<std::array<REAL, 9>,9> M1d({{
            {0  , 0  , 0 , 0  , 0  , 0  ,L+2*M, 0  , 0},
            {0  , 0  , 0 , 0  , 0  , 0  ,L, 0  , 0},
            {0  , 0  , 0 , 0  , 0  , 0  ,L, 0  , 0},
            {0  , 0  , 0 , 0  , 0  , 0  ,0, M  , 0},
            {0  , 0  , 0 , 0  , 0  , 0  ,0, 0  , M},
            {0  , 0  , 0 , 0  , 0  , 0  ,0, 0  , 0},
            {REAL(1)/r  , 0  , 0 , 0  , 0  , 0  ,0, 0  , 0},
            {0  , 0  , 0 , REAL(1)/r  , 0  , 0  ,0, 0  , 0},
            {0  , 0  , 0 , 0, REAL(1)/r  , 0  , 0  ,0,  0}
          }});
    MVFunction<3,9,9,REAL> M1(R2M<3,9,9,REAL>(M1d));
    std::array<std::array<REAL, 9>,9> M2d({{
            {0  , 0  , 0 , 0  , 0  , 0  ,0, L  , 0},
            {0  , 0  , 0 , 0  , 0  , 0  ,0, L+2*M  , 0},
            {0  , 0  , 0 , 0  , 0  , 0  ,0, L  , 0},
            {0  , 0  , 0 , 0  , 0  , 0  ,M, 0  , 0},
            {0  , 0  , 0 , 0  , 0  , 0  ,0, 0  , 0},
            {0  , 0  , 0 , 0  , 0  , 0  ,0, 0  , M},
            {0  , 0  , 0 , REAL(1)/r  , 0  , 0  ,0, 0  , 0},
            {0  , REAL(1)/r  , 0 , 0  , 0  , 0  ,0, 0  , 0},
            {0  , 0  , 0 , 0,0, REAL(1)/r   , 0  ,0,  0}
          }});
    MVFunction<3,9,9,REAL> M2(R2M<3,9,9,REAL>(M2d));
   
    std::array<std::array<REAL, 9>,9> M3d({{
            {0  , 0  , 0 , 0  , 0  , 0  ,0, 0  , L+2*M},
            {0  , 0  , 0 , 0  , 0  , 0  ,0, 0  , L},
            {0  , 0  , 0 , 0  , 0  , 0  ,0, 0  , L},
            {0  , 0  , 0 , 0  , 0  , 0  ,0, 0  , 0},
            {0  , 0  , 0 , 0  , 0  , 0  ,0, 0  , 0},
            {0  , 0  , 0 , 0  , 0  , 0  ,M, 0  , 0},
            {0  , 0  , 0 ,0, REAL(1)/r  ,  0  ,0, 0  , 0},
            {0  , 0  , 0 , 0  , 0  , REAL(1)/r  ,0, 0  , 0},
            {0  , 0  , REAL(1)/r , 0,0, 0   , 0  ,0,  0}
          }});
    MVFunction<3,9,9,REAL> M3(R2M<3,9,9,REAL>(M3d));
   auto d1 = DifferentialOperator<3,9,REAL>({1,0,0},M1);
   auto d2 = DifferentialOperator<3,9,REAL>({0,1,0},M2);
   auto d3 = DifferentialOperator<3,9,REAL>({0,0,1},M3);
   return d1+d2+d3;
  }

  ConstantDifferentialOperator<3,9,REAL> linear_elasticity_const(const REAL& r, const REAL& L, const REAL& M){
    auto M1 = Matrix<9,9,REAL>({{
            {0  , 0  , 0 , 0  , 0  , 0  ,L+2*M, 0  , 0},
            {0  , 0  , 0 , 0  , 0  , 0  ,L, 0  , 0},
            {0  , 0  , 0 , 0  , 0  , 0  ,L, 0  , 0},
            {0  , 0  , 0 , 0  , 0  , 0  ,0, M  , 0},
            {0  , 0  , 0 , 0  , 0  , 0  ,0, 0  , M},
            {0  , 0  , 0 , 0  , 0  , 0  ,0, 0  , 0},
            {REAL(1)/r  , 0  , 0 , 0  , 0  , 0  ,0, 0  , 0},
            {0  , 0  , 0 , REAL(1)/r  , 0  , 0  ,0, 0  , 0},
            {0  , 0  , 0 , 0, REAL(1)/r  , 0  , 0  ,0,  0}
          }});
    auto M2= Matrix<9,9,REAL>({{
            {0  , 0  , 0 , 0  , 0  , 0  ,0, L  , 0},
            {0  , 0  , 0 , 0  , 0  , 0  ,0, L+2*M  , 0},
            {0  , 0  , 0 , 0  , 0  , 0  ,0, L  , 0},
            {0  , 0  , 0 , 0  , 0  , 0  ,M, 0  , 0},
            {0  , 0  , 0 , 0  , 0  , 0  ,0, 0  , 0},
            {0  , 0  , 0 , 0  , 0  , 0  ,0, 0  , M},
            {0  , 0  , 0 , REAL(1)/r  , 0  , 0  ,0, 0  , 0},
            {0  , REAL(1)/r  , 0 , 0  , 0  , 0  ,0, 0  , 0},
            {0  , 0  , 0 , 0,0, REAL(1)/r   , 0  ,0,  0}
          }});
   
    auto M3= Matrix<9,9,REAL>({{
            {0  , 0  , 0 , 0  , 0  , 0  ,0, 0  , L+2*M},
            {0  , 0  , 0 , 0  , 0  , 0  ,0, 0  , L},
            {0  , 0  , 0 , 0  , 0  , 0  ,0, 0  , L},
            {0  , 0  , 0 , 0  , 0  , 0  ,0, 0  , 0},
            {0  , 0  , 0 , 0  , 0  , 0  ,0, 0  , 0},
            {0  , 0  , 0 , 0  , 0  , 0  ,M, 0  , 0},
            {0  , 0  , 0 ,0, REAL(1)/r  ,  0  ,0, 0  , 0},
            {0  , 0  , 0 , 0  , 0  , REAL(1)/r  ,0, 0  , 0},
            {0  , 0  , REAL(1)/r , 0,0, 0   , 0  ,0,  0}
          }});
   auto d1 = ConstantDifferentialOperator<3,9,REAL>({1,0,0},M1);
   auto d2 = ConstantDifferentialOperator<3,9,REAL>({0,1,0},M2);
   auto d3 = ConstantDifferentialOperator<3,9,REAL>({0,0,1},M3);
   return d1+d2+d3;
  }

  DifferentialOperator<3,9,REAL> linear_elasticity_trig(){
    std::shared_ptr<Cinfinity<3,REAL>> z(new Polynomial<3,REAL>(0));
    auto s2x = sine<3>(0,2);
    auto s3y = sine<3>(1,3);
    auto s4z = sine<3>(1,4);
    auto sx = sine<3>(0,1);
    auto L = s2x*s3y*s4z;
    auto M = s2x+s3y;
    auto r = s4z*(sx+s3y);
    MVFunction<3,9,9,REAL> M1({{
            {z  , z  , z , z , z  , z  , L+M+M,   z  , z},
            {z  , z  , z , z , z  , z  , L,     z  , z},
            {z  , z  , z , z , z  , z  , L,     z  , z},
            {z  , z  , z , z , z  , z  , z,     M  , z},
            {z  , z  , z , z , z  , z  , z,     z  , M},
            {z  , z  , z , z , z  , z  , z,     z  , z},
            {r  , z  , z , z , z  , z  , z,     z  , z},
            {z  , z  , z , r , z  , z  , z,     z  , z},
            {z  , z  , z , z , r  , z  , z,     z  , z}
          }});
    MVFunction<3,9,9,REAL> M2({{
            {z  , z  , z , z  , z  , z  ,z, L  , z},
            {z  , z  , z , z  , z  , z  ,z, L+M+M  , z},
            {z  , z  , z , z  , z  , z  ,z, L  , z},
            {z  , z  , z , z  , z  , z  ,M, z  , z},
            {z  , z  , z , z  , z  , z  ,z, z  , z},
            {z  , z  , z , z  , z  , z  ,z, z  , M},
            {z  , z  , z , r  , z  , z  ,z, z  , z},
            {z  , r  , z , z  , z  , z  ,z, z  , z},
            {z  , z  , z , z  , z  , r  ,z, z  ,  z}
          }});
   
    MVFunction<3,9,9,REAL> M3({{
            {z  , z  , z , z  , z  , z  ,z, z  , L+M+M},
            {z  , z  , z , z  , z  , z  ,z, z  , L},
            {z  , z  , z , z  , z  , z  ,z, z  , L},
            {z  , z  , z , z  , z  , z  ,z, z  , z},
            {z  , z  , z , z  , z  , z  ,z, z  , z},
            {z  , z  , z , z  , z  , z  ,M, z  , z},
            {z  , z  , z , z  , r  , z  ,z, z  , z},
            {z  , z  , z , z  , z  , r  ,z, z  , z},
            {z  , z  , r , z  , z  , z  ,z ,z  , z}
          }});
   auto d1 = DifferentialOperator<3,9,REAL>({1,0,0},M1);
   auto d2 = DifferentialOperator<3,9,REAL>({0,1,0},M2);
   auto d3 = DifferentialOperator<3,9,REAL>({0,0,1},M3);
   return d1+d2+d3;
  }
}
#endif
