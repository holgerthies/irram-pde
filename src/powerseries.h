#ifndef POWERSERIES_H
#define POWERSERIES_H
#include <iRRAM.h>
#include "cinfinity.h"
namespace iRRAM{

  REAL error_to_real(const REAL& x)
  {
    sizetype r;
    x.geterror(r);
    REAL err = scale(REAL(signed(r.mantissa)),r.exponent);
    sizetype_exact(err.error);
    return err;
  }
  sizetype real_to_error(const REAL& error)
  {
    sizetype ans;
    REAL zero_one=0;
    sizetype l;
    sizetype_set(l,1,0);
    zero_one.seterror(l);
    REAL e = zero_one*error;
    e.geterror(ans);
    return ans;
    
  }

  class FactorCache{
  private:
    FactorCache() = default;
    ~FactorCache() = default;
    static FactorCache* instance;
    std::vector<std::vector<REAL>> factor_cache;
  public:
    static void init(){
      instance = new FactorCache();
    }

    REAL factor(int n, int m){
      while(factor_cache.size() <= n){
        factor_cache.push_back({1});
      }
      while(factor_cache[n].size() <= m){
        int j = factor_cache[n].size();
        factor_cache[n].push_back(factor_cache[n][j-1]*REAL(n-j+1));
      }
      return factor_cache[n][m];
    }
    static REAL get_factor(int n, int m){
      return instance->factor(n,m);
    }
    static auto get_size(){
      return instance->factor_cache.size();
    }
  };
  FactorCache* FactorCache::instance = 0;
  template<unsigned int d, class T>
  class Powerseries : public Cinfinity<d,T> {
  private:
    std::array<T,1> center;
    REAL radius;
    std::function<T(unsigned int)> coeff_fun;
    std::function<REAL(unsigned int)> bound_fun;
    mutable std::vector<T> coeffs;
    mutable std::vector<REAL> bound;
    REAL get_factor(int n, int m) const{
      return FactorCache::get_factor(n,m);
    }
  protected:
    virtual T get_coefficient_raw(const unsigned int i) const{
      return coeff_fun(i);
    }
    virtual REAL get_bound_raw(const unsigned int i) const{
      return bound_fun(i);
    }
    unsigned int cache_size() const{
      return coeffs.size();
    }
  public:
    Powerseries<d,T>(const std::shared_ptr<Cinfinity<d,T>>& f, const std::array<T,1>& center, const REAL& radius) : center(center), radius(radius) {
      //this->f = [this] (auto ind, auto x){ return sum(ind,x);};
      //this->bound = [this] (auto ind, auto x, auto eps){return 0;};
      coeff_fun = [f,center] (unsigned int i) {return inv_factorial(i)*f->evaluate({i},center);};
      bound_fun = [f,center,radius] (unsigned int i) {return f->get_bound({i},center,radius);};
      
    }
    Powerseries<d,T>() : Powerseries(std::make_shared<Cinfinity<d,T>>(), {0},0){}

    Powerseries<d,T>(const std::function<T(unsigned int)>& coeff_fun, const std::function<REAL(unsigned int)>& bound_fun, const std::array<T,1>& center, const REAL& radius) : center(center), radius(radius), coeff_fun(coeff_fun), bound_fun(bound_fun) {}

    T get_coefficient(const unsigned int i) const{
      auto sz=coeffs.size();
      if(i < sz) return coeffs[i];
      int e = max(i,0);
      coeffs.resize(e+1);
      for(int j=sz; j<=e;j++ ){
        coeffs[j] = get_coefficient_raw(j);
      }
      return coeffs[i];
    }

    REAL get_bound(unsigned int i) const{
    auto sz=bound.size();
    if(i < sz) return bound[i];
    int e = max(i,0);
    bound.resize(e+1);
    for(int j=sz; j<=e;j++){
      bound[j] = get_bound_raw(j);
    }
    return bound[i];
  }

  REAL get_bound(const std::array<unsigned int, d>& index, const std::array<T,d>& x, const REAL& eps) const override{
    return get_bound(index[0]);
  }

  std::array<T,1> get_center() const {return center;}
  REAL get_radius() const {return radius;}
  T evaluate(const std::array<unsigned int,d>& index, const std::array<T,d>& x) const override{
    return sum(index,x);
  }
  T sum(const std::array<unsigned int,d>& index, const std::array<T,d>& x) const {
    T xp=1, xc = (x[0] - center[0]);
    int j = index[0];
    T sum=get_coefficient(j)*xp;
    T ans = sum;
    REAL error = get_bound(j+1)*inv_factorial(j+1)*xc;
    sizetype trunc_error = real_to_error(error), sum_error,total_error;
    sum.geterror(sum_error);
    sizetype_add(total_error, sum_error, trunc_error);
    int i=0;
    while (sizetype_less(sum_error, trunc_error) &&
           (trunc_error.exponent >= ACTUAL_STACK.actual_prec) ){
      i++;
      xp *= xc;
      sum += get_factor(i+j,j)*get_coefficient(i+j)*xp;
      error = get_bound(i+j+1)*inv_factorial(i+j+1)*xp*xc;
      trunc_error = real_to_error(error);
      sum.geterror(sum_error);
      sizetype curr_error;
      sizetype_add(curr_error, sum_error, trunc_error);
      if(sizetype_less(curr_error, total_error)){
        ans = sum;
        ans.seterror(curr_error);
        total_error = curr_error;
      }
    }
    return ans;
  }

  Powerseries<d,T>& operator+=(const Powerseries<d,T>& f){
    *this = *this + f;
    return *this;
  };
};

template<unsigned int d, class T>
class PS_ptr {
private :
  std::shared_ptr<Powerseries<d,T>> ps;
public:  
  PS_ptr() : ps(std::make_shared<Powerseries<d,T>>()) {};
  PS_ptr(const std::shared_ptr<Powerseries<d,T>>& ps) : ps(ps) {};
  Powerseries<d,T>* operator->() const{
    return ps.get();
  }
  T operator()(const std::array<T,d>& x) const {
    return (*ps)(x);
  }
  std::shared_ptr<Cinfinity<d,T>> fun() const{
    return ps;
  }
  PS_ptr<d,T>& operator+=(const PS_ptr<d,T>& p2){
    auto sum = *this + p2;
    *this = sum;
    return *this;
  }
};

template<unsigned int d, class T>
PS_ptr<d,T> operator*(const T& lhs, const PS_ptr<d,T>& rhs){
  return std::make_shared<Powerseries<d,T>>(lhs*rhs.fun(), rhs->get_center(),rhs->get_radius());
}
template<unsigned int d, class T>
PS_ptr<d,T> operator*(const PS_ptr<d,T>& lhs, const PS_ptr<d,T>& rhs){
  REAL new_radius = minimum(lhs->get_radius(), rhs->get_radius());
  // for now assume that the center is the same (this should be
  // fixed later)
  auto new_center = lhs->get_center();
  return std::make_shared<Powerseries<d,T>>(lhs.fun()*rhs.fun(), new_center,new_radius);
}

template<unsigned int d, class T>
PS_ptr<d,T> operator+(const PS_ptr<d,T>& lhs, const PS_ptr<d,T>& rhs){
  REAL new_radius = minimum(lhs->get_radius(), rhs->get_radius());
  // for now assume that the center is the same (this should be
  // fixed later)
  auto new_center = lhs->get_center();
  return std::make_shared<Powerseries<d,T>>(lhs.fun()+rhs.fun(), new_center,new_radius);
}
  
template <unsigned int d, unsigned int m, unsigned int n, class T>
class MVPowerseries : public Matrix<m,n,PS_ptr<d,T>> {
public:
  using Matrix<m,n,PS_ptr<d,T>>::Matrix;
  MVPowerseries(const Matrix<m,n,PS_ptr<d,T>>& M) : Matrix<m,n,PS_ptr<d,T>>(M) {};
  using Matrix<m,n,PS_ptr<d,T>>::operator();

  MVPowerseries<d,m,n,T> operator-();
  MVPowerseries<d,m,n,T> partial_derivative(const unsigned int i) const;
  Matrix<m,n,T> operator()(const std::array<T,d>& x) const{
    Matrix<m,n,T> ans;
    for(unsigned int i=0; i < m; i++){
      for(unsigned int j=0; j < n; j++){
        ans(i,j) = (*this)(i,j)(x);
      }
    }
    return ans;
  }
};
// compute partial derivative
    template <unsigned int d, unsigned int m, unsigned int n, class T>
      MVPowerseries<d,m,n,T> MVPowerseries<d,m,n,T>::partial_derivative(const unsigned int i) const{
      std::array<unsigned int, d> index;
      index[i-1] = 1;
      MVPowerseries<d,m,n,T> ans;
      for(unsigned int i=0; i < m; i++){
        for(unsigned int j=0; j < n; j++){
          ans(i,j) = std::make_shared<Powerseries<d,REAL>>(std::make_shared<Cinfinity<d,REAL>>((*this)(i,j)->derive(index)), (*this)(i,j)->get_center(), (*this)(i,j)->get_radius());;
        }
      }
      return ans;
    }

// redefine overloaded operators on matrix type
    template <unsigned int d, unsigned int m, unsigned int n, class T>
      MVPowerseries<d,m,n,T> operator+(const MVPowerseries<d,m,n,T>& lhs, const MVPowerseries<d,m,n,T>& rhs){
      return add(lhs,rhs);
    }

    template <unsigned int d, unsigned int m, unsigned int n, class T>
      MVPowerseries<d,m,n,T> operator-(const MVPowerseries<d,m,n,T>& lhs, const MVPowerseries<d,m,n,T>& rhs){
      return subtract(lhs,rhs);
    }

    template <unsigned int d, unsigned int m, unsigned int n, unsigned int k, class T>
      MVPowerseries<d,m,k,T> operator*(const MVPowerseries<d,m,n,T>& lhs, const MVPowerseries<d,n,k,T>& rhs){
      return multiply(lhs,rhs);
    }

    template <unsigned int d, unsigned int m, unsigned int n, class T>
      MVPowerseries<d,m,n,T> operator*(const T& lhs, const MVPowerseries<d,m,n,T>& rhs){
      MVPowerseries<d,m,n,T> ans;
      for(int i=0; i<m;i++){
        for(int j=0;j<n;j++){
          ans(i,j) = lhs*rhs(i,j);
        }
      }
      return ans;
    }
    template <unsigned int d, unsigned int m, unsigned int n, class T>
      MVPowerseries<d,m,n,T> operator*(const MVPowerseries<d,m,n,T>& lhs, const T& rhs){
      return multiply(lhs,rhs);
    }
  }

#endif
