#ifndef PDE_H
#define PDE_H
#include <iRRAM.h>
#include "cinfinity.h"
namespace iRRAM{
  // euclidean norm
  template<unsigned int d>
  REAL norm2(const Matrix<d,1,REAL>& v){
    REAL ans=0;
    for(int i = 0; i<d; i++){
      ans += v(i,0)*v(i,0);
    }
    return sqrt(ans);
  }

  template<unsigned int e, unsigned int d, class T>
  class Pde_solution_series{
  private:
    std::array<MVPowerseries<d,e,e,T>, d> f;
    std::vector<MVPowerseries<d,e,1,T>> v;
    std::array<T,d> x;
    MVPowerseries<d,e,1,T> F;
    unsigned int coeffs=0;
  public:
    T get_coefficient(unsigned int k, unsigned int i) {
      int sz = v.size();
      if(sz <= i) v.resize(i+1);
      for(int j=sz; j<=i; j++){
        for(int k=0; k<d; k++){
          v[j] = v[j] + (1/REAL(j))*f[k]*(v[j-1].partial_derivative(k+1));
        }
      }
      return v[i](x)(k,0);
    }

    REAL get_bound(unsigned int k, unsigned int i) const {
      return 1;
    }

    Pde_solution_series(const std::array<MVPowerseries<d,e,e,REAL>, d>& f, const MVPowerseries<d,e,1,REAL>& v, const std::array<T,d>& x) : f(f), v({v}), x(x) {}
  };

  template<unsigned int e, unsigned int d, class T>
  class Pde_local_solution :public MVPowerseries<1,e,1,T> {
  public: 
    Pde_local_solution(const std::array<MVPowerseries<d,e,e,REAL>, d>& f, const MVPowerseries<d,e,1,REAL>& v, const std::array<T,d>& x) {
      std::array<T,1> center ={{0}};;
      auto radius = (v.get_radius() - norm2(v(x)))/REAL(int(d));
      auto series = std::make_shared<Pde_solution_series<e,d,T>>(f,v,x);
      for(int k=0; k<e;k++){
        (*this)(k,0) = std::make_shared<Powerseries<1,REAL>>(
          [k, series] (const std::array<unsigned int, 1>& index) {return series->get_coefficient(k,index[0]);},
          [k, series] (const std::array<unsigned int, 1>& index) {return series->get_bound(k,index[0]);},
          center, radius
          );
      }
    }
  };
}
#endif
