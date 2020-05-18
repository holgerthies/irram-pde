#ifndef CINFIITY_H
#define CINFIITY_H
#include <functional>
#include <iRRAM.h>
#include <combinatorics.h>
#include <array>
namespace iRRAM{
  template<unsigned int n, class T>
  class Cinfinity{
    using fun_type = std::function<T(const std::array<unsigned int, n>&, const std::array<T, n>&)>;
    using bound_type = std::function<REAL(const std::array<unsigned int, n>&, const std::array<T,n>&, const REAL&)>;
  private:
    const fun_type f;
    const bound_type bound;
  public:
    Cinfinity(const fun_type &f, const bound_type &bound) : f(f), bound(bound) {}
    T operator()(const std::array<T,n>&) const;
    T evaluate(const std::array<unsigned int,n>& index, const std::array<T,n>& x) const{
      return f(index,x);
    }

    REAL get_bound(const std::array<unsigned int, n>& index, const std::array<T,n>& x, const REAL& eps) const{
      return bound(index, x, eps);
    }

    Cinfinity<n,T> operator-();
    Cinfinity<n,T> derive(const std::array<unsigned int, n>& index) const;
    template<unsigned int m,class D>
    friend Cinfinity<m,D> operator+(const Cinfinity<m,D>& lhs, const Cinfinity<m,D>& rhs);
    template<unsigned int m,class D>
    friend Cinfinity<m,D> operator-(const Cinfinity<m,D>& lhs, const Cinfinity<m,D>& rhs);
    template<unsigned int m,class D>
    friend Cinfinity<m,D> operator*(const D& lhs, const Cinfinity<m,D>& rhs);
  };

  // evaluation
  template <unsigned int n, class T>
  T Cinfinity<n,T>::operator()(const std::array<T,n>& x) const{
    return f(std::array<unsigned int, n>{0}, x);
  }

  // derivative
  template <unsigned int n, class T>
  Cinfinity<n,T> Cinfinity<n,T>::derive(const std::array<unsigned int, n>& index) const{
    auto add_index = [index] (const std::array<unsigned int,n>& ind) {
                       std::array<unsigned int,n> new_index;
                       std::transform(index.begin(), index.end(), ind.begin(),new_index.begin(), std::plus<unsigned int>());
                       return new_index;
                     };
    return Cinfinity<n,T>([add_index,this] (auto ind, auto x) {return f(add_index(ind), x);}, [add_index,this] (auto ind, auto x, auto eps) {return bound(add_index(ind), x, eps);});
  }

  // unary -
  template <unsigned int n, class T>
  Cinfinity<n,T> Cinfinity<n,T>::operator-(){
    return Cinfinity<n,T>([this] (auto index, auto x) {return -f(index, x);}, bound);
  }

  // addition
  template <unsigned int n, class T>
  Cinfinity<n,T> operator+(Cinfinity<n,T> const& lhs, Cinfinity<n,T> const& rhs){
    return Cinfinity<n,T>([lhs,rhs] (auto index, auto x) {return lhs.f(index,x)+rhs.f(index,x);}, [lhs,rhs](auto index, auto x, auto eps) {return lhs.bound(index,x,eps)+rhs.bound(index,x,eps);} );
  }

  // subtraction
  template <unsigned int n, class T>
  Cinfinity<n,T> operator-(Cinfinity<n,T> const& lhs, Cinfinity<n,T> const& rhs){
    return Cinfinity<n,T>([lhs,rhs] (auto index, auto x) {return lhs.f(index,x)-rhs.f(index,x);}, [lhs,rhs](auto index, auto x, auto eps) {return lhs.bound(index,x,eps)+rhs.bound(index,x,eps);} );
  }
  // scalar multiplication
  template <unsigned int n, class T>
  Cinfinity<n,T> operator*(const T& lhs, Cinfinity<n,T> const& rhs){
    return Cinfinity<n,T>([lhs,rhs] (auto index, auto x) {return lhs*rhs.f(index,x);}, [lhs,rhs](auto index, auto x, auto eps) {return lhs*rhs.bound(index,x,eps);} );
  }
  template <unsigned int n, class T>
  Cinfinity<n,T> operator*(Cinfinity<n,T> const& lhs, const T& rhs){
    return rhs*lhs;
  }
  // multiplication
  // a single summand of the leibniz rule
  template <unsigned int n, class T>
  T product_rule_term(const Cinfinity<n,T>& f, const Cinfinity<n,T>& g, const std::array<unsigned int,n>& alpha, const std::array<unsigned int,n>& beta, const std::array<T, n>& x){
    std::array<unsigned int,n> betap;
    std::transform(alpha.begin(), alpha.end(), beta.begin(),betap.begin(), std::minus<unsigned int>());
    return choose<n>(alpha,beta)*(f.evaluate(beta, x))*(g.evaluate(betap,x));
  }

  // a bound for a single term in the leibniz rule
  template <unsigned int n, class T>
  T product_rule_term_bound(const Cinfinity<n,T>& f, const Cinfinity<n,T>& g, const std::array<unsigned int,n>& alpha, const std::array<unsigned int,n>& beta, const std::array<T, n>& x, const REAL& eps){
    std::array<unsigned int,n> betap;
    std::transform(alpha.begin(), alpha.end(), beta.begin(),betap.begin(), std::minus<unsigned int>());
    return choose<n>(alpha,beta)*(f.get_bound(beta, x, eps))*(g.get_bound(betap,x, eps));
  }

  // leibniz rule

  template <unsigned int n, class T>
  T product_rule(const Cinfinity<n,T>& f, const Cinfinity<n,T>& g, const std::array<unsigned int,n>& alpha, const std::array<T, n>& x){
    T ans;
    for(auto beta : bounded_count<n>(alpha)){
      ans += product_rule_term<n>(f,g,alpha,beta,x);
    }
    return ans;
  }

  template <unsigned int n, class T>
  T product_rule_bound(const Cinfinity<n,T>& f, const Cinfinity<n,T>& g, const std::array<unsigned int,n>& alpha, const std::array<T, n>& x, const REAL& eps){
    REAL ans;
    for(auto beta : bounded_count<n>(alpha)){
      ans += product_rule_term_bound<n>(f,g,alpha,beta,x,eps);
    }
    return ans;
  }

  template <unsigned int n, class T>
  Cinfinity<n,T> operator*(Cinfinity<n,T> const& lhs, Cinfinity<n,T> const& rhs){
    return Cinfinity<n,T>([lhs,rhs] (auto index, auto x) {return product_rule<n>(lhs, rhs, index, x);}, [lhs,rhs](auto index, auto x, auto eps) {return product_rule_bound<n>(lhs,rhs,index,x,eps);} );
  }
}
#endif
