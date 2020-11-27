#ifndef PDE_EXAMPLES_H
#define PDE_EXAMPLES_H
#include <iRRAM.h>
#include "cinfinity.h"
#include "operator.h"
#include "funops.h"
#include "powerseries.h"
#include "mfun.h"
#include "polynomial.h"
#include "diffop.h"
#include "pde.h"
namespace iRRAM{
  // linear elasticity equation
  DifferentialOperator<3,9,REAL> linear_elasticity(const REAL& r, const REAL& L, const REAL& M){
    MVFunction<3,9,9,REAL> M1(R2M<3,9,9,REAL>({{
            {0  , 0  , 0 , 0  , 0  , 0  ,L+2*M, 0  , 0},
            {0  , 0  , 0 , 0  , 0  , 0  ,L, 0  , 0},
            {0  , 0  , 0 , 0  , 0  , 0  ,L, 0  , 0},
            {0  , 0  , 0 , 0  , 0  , 0  ,0, M  , 0},
            {0  , 0  , 0 , 0  , 0  , 0  ,0, 0  , M},
            {0  , 0  , 0 , 0  , 0  , 0  ,0, 0  , 0},
            {REAL(1)/r  , 0  , 0 , 0  , 0  , 0  ,0, 0  , 0},
            {0  , 0  , 0 , REAL(1)/r  , 0  , 0  ,0, 0  , 0},
            {0  , 0  , 0 , 0, REAL(1)/r  , 0  , 0  ,0,  0}
          }}));
    MVFunction<3,9,9,REAL> M2(R2M<3,9,9,REAL>({{
            {0  , 0  , 0 , 0  , 0  , 0  ,0, L  , 0},
            {0  , 0  , 0 , 0  , 0  , 0  ,0, L+2*M  , 0},
            {0  , 0  , 0 , 0  , 0  , 0  ,0, L  , 0},
            {0  , 0  , 0 , 0  , 0  , 0  ,M, 0  , 0},
            {0  , 0  , 0 , 0  , 0  , 0  ,0, 0  , 0},
            {0  , 0  , 0 , 0  , 0  , 0  ,0, 0  , M},
            {0  , 0  , 0 , REAL(1)/r  , 0  , 0  ,0, 0  , 0},
            {0  , REAL(1)/r  , 0 , 0  , 0  , 0  ,0, 0  , 0},
            {0  , 0  , 0 , 0,0, REAL(1)/r   , 0  ,0,  0}
          }}));
   
    MVFunction<3,9,9,REAL> M3(R2M<3,9,9,REAL>({{
            {0  , 0  , 0 , 0  , 0  , 0  ,0, 0  , L+2*M},
            {0  , 0  , 0 , 0  , 0  , 0  ,0, 0  , L},
            {0  , 0  , 0 , 0  , 0  , 0  ,0, 0  , L},
            {0  , 0  , 0 , 0  , 0  , 0  ,0, 0  , 0},
            {0  , 0  , 0 , 0  , 0  , 0  ,0, 0  , 0},
            {0  , 0  , 0 , 0  , 0  , 0  ,M, 0  , 0},
            {0  , 0  , 0 ,0, REAL(1)/r  ,  0  ,0, 0  , 0},
            {0  , 0  , 0 , 0  , 0  , REAL(1)/r  ,0, 0  , 0},
            {0  , 0  , REAL(1)/r , 0,0, 0   , 0  ,0,  0}
          }}));
   auto d1 = DifferentialOperator<3,9,REAL>({1,0,0},M1);
   auto d2 = DifferentialOperator<3,9,REAL>({0,1,0},M2);
   auto d3 = DifferentialOperator<3,9,REAL>({0,0,1},M3);
   return d1+d2+d3;
  }

  template<unsigned int d>
  DifferentialOperator<d,1,REAL> heat_equation(const REAL& alpha){
    DifferentialOperator<d,1, REAL> ans;
    for(int i=0; i<d; i++){
      Multiindex<d> idx;
      idx[i] = 2;
      DifferentialOperator<d,1,REAL> coeff(idx, R2M<d,1,1,REAL>({{{alpha}}}));
      ans = ans+coeff;
    }
    return ans;
  }
}
#endif
