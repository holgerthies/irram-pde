#ifndef POWERSERIES_H
#define POWERSERIES_H
#include <iRRAM.h>
#include "cinfinity.h"
#include <numeric>
namespace iRRAM{

  // euclidean norm
  template<unsigned int d>
  REAL norm2(const std::array<REAL,d>& x){
    REAL ans=0;
    for(int i = 0; i<d; i++){
      ans += x[i]*x[i];
    }
    return sqrt(ans);
  }
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

  // caches n*(n-1)*..*(n-j+1)
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
  class CoefficientCache {
  private:
    std::vector<CoefficientCache<d-1,T>> cache;
  public:
    T get(typename std::array<unsigned int,d>::const_iterator start){
      return cache[*start].get(start+1);
    }

    T get(const std::array<unsigned int,d>& index){
      return get(index.cbegin());
    }

    void put(typename std::array<unsigned int,d>::const_iterator start, const T& val){
      if(cache.size() <= (*start)+1)
        cache.resize((*start)+1);
      cache[*start].put(start+1,val);
    }

    void put(const std::array<unsigned int,d>& index, const T& val){
      put(index.cbegin(), val);
    }

    bool isvalid(typename std::array<unsigned int,d>::const_iterator start){
      if(cache.size() <= *start) return false;
      return cache[*start].isvalid(start+1);
    }

    bool isvalid(const std::array<unsigned int,d>& index){
      return isvalid(index.cbegin());
    }
  };

  template<class T>
  class CoefficientCache<1,T> {
  private:
    std::vector<T> cache;
  public:
    T get(typename std::array<unsigned int,1>::const_iterator start){
      return cache[*start];
    }

    T get(const std::array<unsigned int,1>& index){
      return get(index.cbegin());
    }

    void put(typename std::array<unsigned int,1>::const_iterator start, const T& val){
      cache.resize(*start+1);
      cache[*start] = val;
    }

    void put(const std::array<unsigned int,1>& index, const T& val){
      put(index.cbegin(), val);
    }

    bool isvalid(typename std::array<unsigned int,1>::const_iterator start){
      return (cache.size() > *start);
    }

    bool isvalid(const std::array<unsigned int,1>& index){
      return isvalid(index.cbegin());
    }
  };

  template<unsigned int d, class T>
  class Powerseries : public Cinfinity<d,T> {
  private:
    std::function<T(const std::array<unsigned int, d>&)> coeff_fun;
    std::function<REAL(const std::array<unsigned int,d>&)> bound_fun;
    mutable std::vector<T> coeffs = {};
    mutable CoefficientCache<d,T> coeff_cache = {};
    mutable CoefficientCache<d,REAL> bound_cache = {};
    REAL get_factor(int n, int m) const{
      return FactorCache::get_factor(n,m);
    }
    REAL get_factor(const std::array<unsigned int,d>& n,const std::array<unsigned int,d>& m) const{
      REAL ans=1;
      for(int i=0; i<d; i++)
        ans *= FactorCache::get_factor(n[i], m[i]);
      return ans;
    }
    REAL get_bound_complex(const std::array<unsigned int,d>& index, const std::array<T,d>& x) const{
      REAL r = get_radius();
      REAL dist = r-norm2<d>(x);
      if(dist < 0){
        cout << "Warning: radius too small" << std::endl;
      }
      unsigned int k = std::accumulate(index.begin(), index.end(), 0);
      return M*power(r, d)/power(dist,k+d);
    }
  protected:
    std::array<T,d> center;
    REAL radius,M;
    virtual T get_coefficient_raw(const std::array<unsigned int, d>& index) const{
      return coeff_fun(index);
    }
    virtual REAL get_bound_raw(const std::array<unsigned int, d>& index) const{
      return bound_fun(index);
    }
    unsigned int cache_size() const{
      return coeffs.size();
    }
  public:
    Powerseries<d,T>(const std::shared_ptr<Cinfinity<d,T>>& f, const std::array<T,d>& center, const REAL& radius, const REAL& M) : center(center), radius(radius), M(M) {
      coeff_fun = [f,center] (const std::array<unsigned int,d>& index) {return inv_factorial<d>(index)*f->evaluate(index,center);};
      bound_fun = [f,center,radius] (const std::array<unsigned int,d>& index) {return f->get_bound({index},center,radius);};
    }
    Powerseries<d,T>() : Powerseries(std::make_shared<Cinfinity<d,T>>(), {0},100000, 0){}
    Powerseries<d,T>(const Powerseries<d,T>& ps) : center(ps.get_center()), radius(ps.get_radius()), M(ps.get_M()){
      coeff_fun = [&ps] (const std::array<unsigned int,d>& index) {return ps.get_coefficient(index);};
      bound_fun = [&ps] (const std::array<unsigned int,d>& index) {return ps.get_bound(index);};
    }

    Powerseries<d,T>(const std::function<T(const std::array<unsigned int,d>&)>& coeff_fun, const std::function<REAL(const std::array<unsigned int, d>&)>& bound_fun, const std::array<T,1>& center, const REAL& radius, const REAL& M) :  coeff_fun(coeff_fun), bound_fun(bound_fun),center(center), radius(radius), M(M) {
    }

    T get_coefficient(const std::array<unsigned int,d>& index) const{
      if(!coeff_cache.isvalid(index)){
        // make sure previous coefficients are already cached
        std::array<unsigned int, d> pindex(index);
        for(int i=0; i<d; i++){
          if(pindex[i] > 0){
            pindex[i]--;
            get_coefficient(pindex);
            pindex[i]++;
          }
        }
        if(!coeff_cache.isvalid(index))
          coeff_cache.put(index,get_coefficient_raw(index));
      }
      return coeff_cache.get(index);
    }

    REAL get_bound(const std::array<unsigned int,d>& index) const{
      if(!bound_cache.isvalid(index)){
        // make sure previous coefficients are already cached
        std::array<unsigned int, d> pindex(index);
        for(int i=0; i<d;i++){
          if(pindex[i] > 0){
            pindex[i]--;
            get_bound(pindex);
            pindex[i]++;
          }
        }
        bound_cache.put(index,get_bound_raw(index));
      }
      return bound_cache.get(index);
    }

    REAL get_bound(const std::array<unsigned int, d>& index, const std::array<T,d>& x, const REAL& eps) const override{
      return get_bound(index);
    }

    std::array<T,d> get_center() const {return center;}
    REAL get_radius() const {return radius;}
    REAL get_M() const {return M;}
    T evaluate(const std::array<unsigned int,d>& index, const std::array<T,d>& x) const override{
      return sum(index,x);
    }

    T sum(const std::array<unsigned int,d>& index, const std::array<T,d>& x) const {
      std::array<T, d> xc;
      std::transform(x.begin(), x.end(), center.begin(),xc.begin(), std::minus<T>());
      T sum=get_factor(index,index)*get_coefficient(index);
      T ans = sum;
      REAL error = 0;
      for(auto p : partitions<d>(1)){
        std::array<unsigned int,d> p2;
        std::transform(index.begin(), index.end(), p.begin(),p2.begin(), std::plus<unsigned int>());
        REAL bnd = minimum(get_bound(p2), get_bound_complex(p2, x));
        error += bnd*power<d,T>(xc,p2);
      }
      int i = 0;
      sizetype trunc_error = real_to_error(error), sum_error,total_error;
      sum.geterror(sum_error);
      sizetype_add(total_error, sum_error, trunc_error);
      while (sizetype_less(sum_error, trunc_error) &&
             (trunc_error.exponent >= ACTUAL_STACK.actual_prec) ){
        i++;
        for(auto p : partitions<d>(i)){
          std::array<unsigned int,d> p2;
          std::transform(index.begin(), index.end(), p.begin(),p2.begin(), std::plus<unsigned int>());
          sum += get_factor(p2,index)*get_coefficient(p2)*power<d,T>(xc, p2);
        }
        error = 0;
        for(auto p : partitions<d>(i+1)){
          std::array<unsigned int,d> p2;
          std::transform(index.begin(), index.end(), p.begin(),p2.begin(), std::plus<unsigned int>());
          REAL bnd = minimum(inv_factorial<d>(p)*get_bound(p2), get_bound_complex(p2, x));
          error += bnd*power<d,T>(xc,p2);
        }
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
    }

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
  CinfinityPtr<d,T> toPowerseries(const CinfinityPtr<d,T>& f, const std::array<T,d> center, const REAL& radius, const REAL& M){
    return std::make_shared<Powerseries<d,T>>(f, center, radius,M);
  }

  template<unsigned int d, class T>
  PS_ptr<d,T> operator*(const T& lhs, const PS_ptr<d,T>& rhs){
    return std::make_shared<Powerseries<d,T>>(lhs*rhs.fun(), rhs->get_center(),rhs->get_radius(), rhs->get_M()*lhs);
  }
  template<unsigned int d, class T>
PS_ptr<d,T> operator*(const PS_ptr<d,T>& lhs, const PS_ptr<d,T>& rhs){
  REAL new_radius = minimum(lhs->get_radius(), rhs->get_radius());
  // for now assume that the center is the same (this should be
  // fixed later)
  auto new_center = rhs->get_center();
  return std::make_shared<Powerseries<d,T>>(lhs.fun()*rhs.fun(), new_center,new_radius, rhs->get_M()*lhs->get_M());
}

template<unsigned int d, class T>
PS_ptr<d,T> operator+(const PS_ptr<d,T>& lhs, const PS_ptr<d,T>& rhs){
  REAL new_radius = minimum(lhs->get_radius(), rhs->get_radius());
  // for now assume that the center is the same (this should be
  // fixed later)
  auto new_center = rhs->get_center();
  return std::make_shared<Powerseries<d,T>>(lhs.fun()+rhs.fun(), new_center,new_radius, rhs->get_M()+lhs->get_M());
}
  
template <unsigned int d, unsigned int m, unsigned int n, class T>
class MVPowerseries : public Matrix<m,n,PS_ptr<d,T>> {
private:
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

  auto get_center() const{
    return (*this)(0,0)->get_center();
  }

  REAL get_radius() const{
    REAL ans = (*this)(0,0)->get_radius();
    for(int i=0; i<m;i++){
      for(int j=0; j<n; j++){
        ans = minimum(ans, (*this)(i,j)->get_radius());
      }
    }
    return ans;
  }

  REAL get_bound(const std::array<unsigned int,d>& index) const{
    REAL ans = (*this)(0,0)->get_bound(index);
    for(int i=0; i<m;i++){
      for(int j=0; j<n; j++){
        ans = maximum(ans, (*this)(i,j)->get_bound(index));
      }
    }
    return ans;
  }

};
// compute partial derivative
  template <unsigned int d, unsigned int m, unsigned int n, class T>
  MVPowerseries<d,m,n,T> MVPowerseries<d,m,n,T>::partial_derivative(const unsigned int i) const{
    std::array<unsigned int, d> index{};
    index[i-1] = 1;
    MVPowerseries<d,m,n,T> ans;
    for(unsigned int i=0; i < m; i++){
      for(unsigned int j=0; j < n; j++){
        ans(i,j) = std::make_shared<Powerseries<d,REAL>>(std::make_shared<Cinfinity<d,REAL>>((*this)(i,j)->derive(index)), (*this)(i,j)->get_center(), (*this)(i,j)->get_radius(), (*this)(i,j)->get_M());
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

template<unsigned int d, unsigned int m, unsigned int n,class T>
MVPowerseries<d,m,n,T> toPowerseries(const MVFunction<d,m,n,T>& f, const std::array<T,d> center, const REAL& radius, const REAL& M){
  MVPowerseries<d,m,n,T> ans;
  for(int i=0; i<m;i++){
    for(int j=0;j<n;j++){
      ans(i,j) = std::make_shared<Powerseries<d,T>>(f(i,j), center, radius,M);
    }
  }
  return ans;
}
}


#endif
