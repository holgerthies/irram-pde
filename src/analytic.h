#ifndef ANALYTIC_H
#define ANALYTIC_H
#include "powerseries.h"
#include <numeric>
namespace iRRAM{

  template<unsigned int m, unsigned int d, class T>
  struct PS_vec : public std::vector<PS_vec<m-1,d,T>> {
    PS_vec(int n = 0) : std::vector<PS_vec<m-1,d,T>>(n, PS_vec<m-1,d,T>(n)) {}
  };
  template<unsigned int d, typename T>
  struct PS_vec<1,d, T> : public std::vector<PS_ptr<d,T>> {
    PS_vec(int n = 0) : std::vector<PS_ptr<d,T>>(n) {}
  };

  // Function analytic on [0,1]^d given by covering of powerseries
  template<unsigned int d, class T>
  class ANALYTIC {
  private:
     int L,M;
  public:
    PS_vec<d,d,T> ps;
    ANALYTIC(const int L, const int M) : L(L), M(M), ps(PS_vec<d,d,T>(2*L+1)) {}
    int get_L() const {return L;}
    int get_M() const {return M;}



  };

  template<unsigned int m, unsigned int d, class T>
  struct to_analytic_recursor{
    std::shared_ptr<Cinfinity<d,REAL>> f;
    PS_vec<m,d,T> to_analytic(const int L, const int M, const std::vector<T>& center){
        to_analytic_recursor<m-1,d,T> r;
        r.f = f.deep_copy();
        REAL ci = REAL(1)/REAL(2*L);
        PS_vec<m,d,T> ans;
        for(int i=0; i<L; i++){
            std::vector<T> new_center(center);
            new_center.push_back(ci);
            ans.push_back(r.to_analytic(L, M,new_center));
            ci += REAL(1)/REAL(L);
        } 
        return ans;
    } 
  };

  template<unsigned int d, class T>
  struct to_analytic_recursor<1,d,T>{
    PS_vec<1,d,T> to_analytic(std::shared_ptr<Cinfinity<d,REAL>> f, const int L, const int M, const std::vector<T>& center){
        REAL ci = REAL(1)/REAL(2*L);
        PS_vec<1,d,T> ans(L);
        REAL r = REAL(1)/REAL(L);
        for(int i=0; i<L; i++){
            vector<REAL, d> new_center;
            for(int j=0; j<d-1;j++){
                new_center[j] = center[j]; 
            }
            new_center[d-1] = ci;
            std::shared_ptr<Cinfinity<d,REAL>> g = f->deep_copy();
            g->set_center(new_center);
            auto ps = std::make_shared<Powerseries<d,T>>([g] (const Multiindex<d>& index) {return g->get_derivative(index);}, new_center, r, M);
            ans[i] = PS_ptr<d,T>(ps);
            ci += REAL(1)/REAL(L);
        }
        return ans;
    }
  };

  template<unsigned int d, class T>
  ANALYTIC<d,T> to_analytic(const CinfinityPtr<d,T> f, const int L, const int M){
      to_analytic_recursor<d,d,T> r;
      ANALYTIC<d,T> ans(L,M);
      ans.ps = r.to_analytic(f, L, M, std::vector<T>());
      return ans;
  }
}

#endif
