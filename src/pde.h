#ifndef PDE_H
#define PDE_H
#include "cinfinity.h"
#include "powerseries.h"
#include "diffop.h"
namespace iRRAM{
  template<unsigned int d, unsigned int e,class C, class T>
  class Pde_solution_series{
  private:
    AbstractDifferentialOperator<d,e,C, T> D;
    std::vector<MVFunction<d,e,1,T>> v;
  public:
    T get_coefficient(unsigned int m, unsigned int i) {
      int sz = v.size();
      if(sz <= i) v.resize(i+1);
      for(int j=sz; j<=i; j++){
          v[j] = (1/REAL(j))*D(v[j-1]);
      }
      return v[i](m,0)->get_derivative({0});
    }

    Pde_solution_series(const AbstractDifferentialOperator<d,e,C,T>& D, const MVFunction<d,e,1,REAL>& v, const vector<T,d>& x) : D(D), v({v}) {
      this->v[0].set_center(x);
    }
  };

  template<unsigned int d, unsigned int e,class T>
  class Pde_solution_series_constant{
  private:
    ConstantDifferentialOperator<d,e,T> D;
    MVFunction<d,e,1,REAL> v;
    std::vector<ConstantDifferentialOperator<d,e,T>> Dn;
  public:
    T get_coefficient(unsigned int m, unsigned int i) {
      int sz = Dn.size();
      if(sz <= i) Dn.resize(i+1);
      for(int j=sz; j<=i; j++){
          Dn[j] = D*Dn[j-1];
      }
      return inv_factorial(i)*Dn[i](v)(m,0)->get_derivative({0});
    }

    Pde_solution_series_constant(const ConstantDifferentialOperator<d,e,T>& D, const MVFunction<d,e,1,REAL>& v, const vector<T,d>& x) : D(D), v(v), Dn({D}) {
      this->v.set_center(x);
    }
  };


  REAL get_sol_r(unsigned int d, unsigned int e, const REAL& r, const REAL& M){
    return r/(5*int(d)*int(d+1)*int(e)*M);
  }
  REAL get_sol_r_const(unsigned int d, unsigned int e, const REAL& r, const REAL& M){
    return r/(int(d)*int(e)*M);
  }
  template<unsigned int e, unsigned int d, class T>
  std::array<Powerseries<1,T>, e> solve_pde(const DifferentialOperator<d,e,T>& D, const MVFunction<d,e,1,REAL>& v, const vector<T,d>& x, const REAL& r, const REAL& M) {
    auto series = std::make_shared<Pde_solution_series<d,e,DifferentialOperatorCoefficient<d,e,T>, T>>(D, v, x);
    std::array<Powerseries<1,T>,e> ans;
    REAL r_new = get_sol_r(d, e, r, M);
    for(int i=0; i<e; i++){
      ans[i] = Powerseries<1,T>([i, series] (const Multiindex<1>& index) {return series->get_coefficient(i, index[0]);}, {0}, r_new, M);
    }
    return ans;
  }

  template<unsigned int e, unsigned int d, class T>
  std::array<Powerseries<1,T>, e> solve_pde_constant(const ConstantDifferentialOperator<d,e,T>& D, const MVFunction<d,e,1,REAL>& v, const vector<T,d>& x, const REAL& r, const REAL& M) {
    auto series = std::make_shared<Pde_solution_series<d,e,DifferentialOperatorConstantCoefficient<d,e,T>, T>>(D, v, x);
    std::array<Powerseries<1,T>,e> ans;
    REAL r_new = get_sol_r_const(d, e, r, M);
    for(int i=0; i<e; i++){
      ans[i] = Powerseries<1,T>([i, series] (const Multiindex<1>& index) {return series->get_coefficient(i, index[0]);}, {0}, r_new, M);
    }
    return ans;
  }

  template<unsigned int e, unsigned int d, class T>
  std::array<Powerseries<1,T>, e> solve_pde_constant_method2(const ConstantDifferentialOperator<d,e,T>& D, const MVFunction<d,e,1,REAL>& v, const vector<T,d>& x, const REAL& r, const REAL& M) {
    auto series = std::make_shared<Pde_solution_series_constant<d,e,T>>(D, v, x);
    std::array<Powerseries<1,T>,e> ans;
    REAL r_new = get_sol_r_const(d, e, r, M);
    for(int i=0; i<e; i++){
      ans[i] = Powerseries<1,T>([i, series] (const Multiindex<1>& index) {return series->get_coefficient(i, index[0]);}, {0}, r_new, M);
    }
    return ans;
  }
  // template<unsigned int e, unsigned int d, class T>
  // class Pde_solution_series{
  // private:
  //   std::array<MVPowerseries<d,e,e,T>, d> f;
  //   std::vector<MVPowerseries<d,e,1,T>> v;
  //   std::array<T,d> x;
  //   MVPowerseries<d,e,1,T> F;
  //   unsigned int coeffs=0;
  // public:
  //   T get_coefficient(unsigned int k, unsigned int i) {
  //     int sz = v.size();
  //     if(sz <= i) v.resize(i+1);
  //     for(int j=sz; j<=i; j++){
  //       for(int k=0; k<d; k++){
  //         v[j] = v[j] + (1/REAL(j))*f[k]*(v[j-1].partial_derivative(k+1));
  //       }
  //     }
  //     return v[i](x)(k,0);
  //   }

  //   REAL get_bound(unsigned int k, unsigned int i) const {
  //     return 1;
  //   }

  //   Pde_solution_series(const std::array<MVPowerseries<d,e,e,REAL>, d>& f, const MVPowerseries<d,e,1,REAL>& v, const std::array<T,d>& x) : f(f), v({v}), x(x) {}
  // };

  // template<unsigned int e, unsigned int d, class T>
  // class Pde_local_solution :public MVPowerseries<1,e,1,T> {
  // public: 
  //   Pde_local_solution(const std::array<MVPowerseries<d,e,e,REAL>, d>& f, const MVPowerseries<d,e,1,REAL>& v, const std::array<T,d>& x) {
  //     std::array<T,1> center ={{0}};
  //     std::array<unsigned int,d> z{};
  //     auto fradius = f[0].get_radius() - norm2(v(x));
  //     auto fM = f[0].get_bound(z);
  //     auto r = (fradius/(REAL(int(d))*fM));
  //     auto M = fradius;
  //     auto series = std::make_shared<Pde_solution_series<e,d,T>>(f,v,x);
  //     for(int k=0; k<e;k++){
  //       (*this)(k,0) = std::make_shared<Powerseries<1,REAL>>(
  //         [k, series] (const std::array<unsigned int, 1>& index) {return series->get_coefficient(k,index[0]);},
  //         [r,M] (const std::array<unsigned int, 1>& index) {return M/power(r,index[0]);},
  //         center, r,M
  //         );
  //     }
  //   }
    //   Pde_local_solution(const std::array<MVFunction<d,e,e,T>, d>& f, const MVFunction<d,e,1,T>& v, const std::array<T,d>& x) {
    //     auto vx = v(x);
    //     std::array<T,d> center;
    //     for(int i=0; i<d; i++)
    //       center[i] = vx(i,0);
    //     std::array<MVPowerseries<d,e,e,T>, d> ps;
    //     MVPowerseries<d,e,1,T> vps = toPowerseries<d,e,1,T>(v,center,1);
    //     for(int i=0; i<d; i++)
    //       ps[i] = toPowerseries<d,e,e,T>(f[i], center,1);
    //     auto fradius = 1;
    //     std::array<unsigned int,d> z{};
    //     auto fM = maximum(ps[0].get_bound(z),1);
    //     auto r = (fradius/(REAL(int(d))*fM));
    //     auto M = fradius;
    //     auto series  = std::make_shared<Pde_solution_series<e,d,T>>(ps,vps,x);
    //     std::array<T,1> zero ={{0}};
    //     for(int k=0; k<e;k++){
    //       (*this)(k,0) = std::make_shared<Powerseries<1,REAL>>(
    //         [k, series] (const std::array<unsigned int, 1>& index) {return series->get_coefficient(k,index[0]);},
    //         [r,M] (const std::array<unsigned int, 1>& index) {return M/power(r,index[0]);},
    //         zero, r
    //         );
    //     }
    //   }
}
#endif
