#include <iRRAM.h>
#include "cinfinity.h"
#include "operator.h"
#include "funops.h"
// #include "matrix.h"
#include "fourier.h"
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

ANALYTIC<1,REAL> heat_equation(const ANALYTIC<1,REAL>& v, const REAL& t){ 
    int L = v.get_L();
    int M = v.get_M();
    Analytic_to_Fourier F(L,v);
    Fourier_to_Analytic G(L,M, [F,t] (const int k) {return F.get_coeff(k)*exp(-4*pi()*pi()*k*k*t); });
    ANALYTIC<1,REAL> u(v.get_L(), v.get_M());
    for(int i=0; i<L; i++){
        vector<REAL,1> center;
        center[0] = 0;
        auto ps = std::make_shared<Powerseries<1,REAL>>([G,i] (const Multiindex<1>& index) {return G.get_ps_coeff(i,index[0]);}, center, L, M);
        u.ps[i] = ps;
    }
    return u;
}


void compute(){
  struct rusage usage;
  struct timeval start, end;
  FactorCache::init();  
  auto f = sine<1>(0,2*pi());
  f->set_center({0});
  auto g = derive(f,{4})+f;
  auto fp = to_analytic(g, 4,1);
  int max_iter;
  static int iteration_counter=0;
  iRRAM::cout << "choose number of iterations" << std::endl;
  iRRAM::cin >>  max_iter;

  iteration_counter++; 
  if(iteration_counter == max_iter) return;

  auto u = heat_equation(fp, REAL("0.001"));
//   Analytic_to_Fourier F(10,fp);
//   for(int l=0; l<2; l++){
//       for(int j=0; j < 3; j++){
//         cout << l << " " << j << " " << fp.ps[l]->get_coefficient({j}) << std::endl;
//       }
//   }

//   Fourier_to_Analytic G(10,1, [F] (const int k) {return F.get_coeff(k); });
//   std::vector<COMPLEX> fc(20);
//   COMPLEX ans(0,0);
//   REAL t = 1;
//   REAL x = REAL(1)/REAL(3);
//   for(int i=-20; i<20; i++){
//       COMPLEX fci = F.get_coeff(i);
//       COMPLEX im(0,1);
//       cout << i << ":" << real(fci) << " " << imag(fci) <<  std::endl;
//       ans = ans + fci*exp(REAL(2)*pi()*i*im*x);
//   }
//   cout << real(ans) << " " << imag(ans) << std::endl;
//    g->set_center({x});
//    cout << g->get_derivative({0}) << std::endl;

   getrusage(RUSAGE_SELF, &usage);
   start = usage.ru_utime;
  
   getrusage(RUSAGE_SELF, &usage);
   start = usage.ru_utime;


   REAL y = u.ps[1]->sum({REAL(1)/REAL(2*fp.get_L())});
   cout<<y<<std::endl;
   getrusage(RUSAGE_SELF, &usage);
   end = usage.ru_utime;
   auto iteration_time =  end.tv_sec-start.tv_sec+double(end.tv_usec-start.tv_usec)/1000000;

   sizetype error;
   y.geterror(error);
   int error_exp_normalized;
   unsigned long mantissa = error.mantissa;
   error_exp_normalized = error.exponent;
   while(mantissa > 1){
     mantissa /= 2;
     error_exp_normalized++;
   }
  
  
   std::cout << " time:" << iteration_time;
   std::cout << " iteration:" << iteration_counter;
   std::cout << " precision:" << actual_stack().actual_prec;
   std::cout << " error_mantissa:"  << error.mantissa;
   std::cout << " error_exponent:" << error.exponent;
   std::cout << " normalized:" << error_exp_normalized;
  
   std::cout << std::endl;
//   std::vector<COMPLEX> Xc(20);
//   COMPLEX ans(0,0);
//   REAL x = REAL(1)/REAL(3);
//   int j=1;
//   int L =2;
//   for(int k=1; k<5000; k++){
//       Xc[k] = COMPLEX(0,0);
//       for(int l = 0; l<L; l++){
//           Xc[k] = Xc[k] + F.getX(l,k,j);
//       }
//       //cout << real(Xc[k]) << " " << imag(Xc[k]) << std::endl;
//       COMPLEX im(0,1);
//       ans = ans + Xc[k]*exp(REAL(-2)*pi()*k*im*x);
//   }
//   cout << real(ans) << " " << imag(ans) << std::endl;
if(REAL(1) < REAL(1)) return;
}

