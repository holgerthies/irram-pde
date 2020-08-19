#ifndef POWERSERIES_H
#define POWERSERIES_H
#include <iRRAM.h>
#include "polynomial.h"
#include <numeric>
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



  template<unsigned int d, class T>
  class Powerseries {
  private:
    std::function<T(const Multiindex<d>&)> coeff_fun;
    mutable Polynomial<d,T> trunc_series;
    REAL get_factor(int n, int m) const{
      return FactorCache::get_factor(n,m);
    }
    REAL get_factor(const std::array<unsigned int,d>& n,const std::array<unsigned int,d>& m) const{
      REAL ans=1;
      for(int i=0; i<d; i++)
        ans *= FactorCache::get_factor(n[i], m[i]);
      return ans;
    }
    REAL get_bound_complex(const int deg, const std::array<T,d>& x) const{
      REAL r = get_radius();
      REAL dist = r-inf_norm(x);
      if(dist < 0){
        cout << "Warning: radius too small" << std::endl;
      }
      return M*power(r, d)/power(dist,deg+d);
    }

    void extend_trunc_series() const{
      for(auto index : max_degree_indices<d>(trunc_series.get_degree()+1)){
        trunc_series.set_coefficient(index,coeff_fun(index));
      }     
    }

  protected:
    vector<T,d> center;
    REAL radius,M;
  public:
    Powerseries<d,T>(const Powerseries<d,T>& ps) : center(ps.get_center()), radius(ps.get_radius()), M(ps.get_M()){
      coeff_fun = [&ps] (const std::array<unsigned int,d>& index) {return ps.get_coefficient(index);};
    }

    Powerseries<d,T>(const std::function<T(const Multiindex<d>&)>& coeff_fun, const vector<T,d>& center, const REAL& radius, const REAL& M) :  coeff_fun(coeff_fun), center(center), radius(radius), M(M) {
    }

    T get_coefficient(const Multiindex<d>& index) const{
      while(trunc_series.get_degree() <= max(index)){
        extend_trunc_series();
      }
      return trunc_series.get_derivative(index);
    }

    vector<T,d> get_center() const {return center;}

    REAL get_radius() const {return radius;}

    REAL get_M() const {return M;}

    // T evaluate(const std::array<unsigned int,d>& index, const std::array<T,d>& x) const override{
    //   return sum(index,x);
    // }

    T partial_sum(const vector<T,d>& x, const int degree) const {
      while(trunc_series.get_degree() <= degree){
        extend_trunc_series();
      }
      return trunc_series(x);
    }

    int get_cache_size() const{
      return trunc_series.get_degree();
    }

    Polynomial<d,T> get_truncated_series(const int degree) const{
      while(trunc_series.get_degree() <= degree){
        extend_trunc_series();
      }
      return trunc_series;
    }

    T sum(const vector<T,d>& x) const {
      vector<T, d> xc = x - center;
      auto xmax = inf_norm(xc);
      REAL error_factor = power(xmax/radius,d);
      REAL error = M*power(1-xmax/radius, d);
      Multiindex<d> zero{};
      REAL sum = get_coefficient(zero);
      REAL ans = sum;
      sizetype trunc_error = real_to_error(error), sum_error,total_error;
      sum.geterror(sum_error);
      sizetype_add(total_error, sum_error, trunc_error);
      int deg=0;
      while (sizetype_less(sum_error, trunc_error) &&
             (trunc_error.exponent >= ACTUAL_STACK.actual_prec)){
        deg++;
        for(auto idx : max_degree_indices<d>(deg)){
          sum += get_coefficient(idx)*power(xc, idx);
        }
        error *= error_factor;
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
    // T sum(const std::array<unsigned int,d>& index, const std::array<T,d>& x) const {
    //   std::array<T, d> xc;
    //   std::transform(x.begin(), x.end(), center.begin(),xc.begin(), std::minus<T>());
      
    //   T sum=get_factor(index,index)*get_coefficient(index);
    //   T ans = sum;
    //   REAL error = 0;
    //   for(auto p : partitions<d>(1)){
    //     std::array<unsigned int,d> p2;
    //     std::transform(index.begin(), index.end(), p.begin(),p2.begin(), std::plus<unsigned int>());
    //     REAL bnd = minimum(get_bound(p2), get_bound_complex(p2, x));
    //     error += bnd*power<d,T>(xc,p2);
    //   }
    //   int i = 0;
    //   sizetype trunc_error = real_to_error(error), sum_error,total_error;
    //   sum.geterror(sum_error);
    //   sizetype_add(total_error, sum_error, trunc_error);
    //   while (sizetype_less(sum_error, trunc_error) &&
    //          (trunc_error.exponent >= ACTUAL_STACK.actual_prec) ){
    //     i++;
    //     for(auto p : partitions<d>(i)){
    //       std::array<unsigned int,d> p2;
    //       std::transform(index.begin(), index.end(), p.begin(),p2.begin(), std::plus<unsigned int>());
    //       sum += get_factor(p2,index)*get_coefficient(p2)*power<d,T>(xc, p2);
    //     }
    //     error = 0;
    //     for(auto p : partitions<d>(i+1)){
    //       std::array<unsigned int,d> p2;
    //       std::transform(index.begin(), index.end(), p.begin(),p2.begin(), std::plus<unsigned int>());
    //       REAL bnd = minimum(inv_factorial<d>(p)*get_bound(p2), get_bound_complex(p2, x));
    //       error += bnd*power<d,T>(xc,p2);
    //     }
    //     trunc_error = real_to_error(error);
    //     sum.geterror(sum_error);
    //     sizetype curr_error;
    //     sizetype_add(curr_error, sum_error, trunc_error);
    //     if(sizetype_less(curr_error, total_error)){
    //       ans = sum;
    //       ans.seterror(curr_error);
    //       total_error = curr_error;
    //     }
    //   }
    //   return ans;
    // }

    // Powerseries<d,T>& operator+=(const Powerseries<d,T>& f){
    //   *this = *this + f;
    //   return *this;
    // }

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
    T operator()(const vector<T,d>& x){
      return ps->sum(x);
    }
    PS_ptr<d,T>& operator+=(const PS_ptr<d,T>& p2){
      auto sum = *this + p2;
      *this = sum;
      return *this;
    }
  };

  template<unsigned int d, class T>
  PS_ptr<d,T> to_powerseries(const CinfinityPtr<d,T>& f, const vector<T,d>& center, const REAL& r, const REAL& M){
    auto g = f->deep_copy();
    g->set_center(center);
    return PS_ptr<d,T>(std::make_shared<Powerseries<d,T>>([g] (const Multiindex<d>& index) {return g->get_derivative(index);}, center, r, M));
  }

// template <unsigned int d, unsigned int m, unsigned int n, class T>
// class MVPowerseries : public Matrix<m,n,PS_ptr<d,T>> {
// private:
// public:
//   using Matrix<m,n,PS_ptr<d,T>>::Matrix;
//   MVPowerseries(const Matrix<m,n,PS_ptr<d,T>>& M) : Matrix<m,n,PS_ptr<d,T>>(M) {};
//   using Matrix<m,n,PS_ptr<d,T>>::operator();

//   MVPowerseries<d,m,n,T> operator-();
//   MVPowerseries<d,m,n,T> partial_derivative(const unsigned int i) const;
//   Matrix<m,n,T> operator()(const std::array<T,d>& x) const{
//     Matrix<m,n,T> ans;
//     for(unsigned int i=0; i < m; i++){
//       for(unsigned int j=0; j < n; j++){
//         ans(i,j) = (*this)(i,j)(x);
//       }
//     }
//     return ans;
//   }

//   auto get_center() const{
//     return (*this)(0,0)->get_center();
//   }

//   REAL get_radius() const{
//     REAL ans = (*this)(0,0)->get_radius();
//     for(int i=0; i<m;i++){
//       for(int j=0; j<n; j++){
//         ans = minimum(ans, (*this)(i,j)->get_radius());
//       }
//     }
//     return ans;
//   }

//   REAL get_bound(const std::array<unsigned int,d>& index) const{
//     REAL ans = (*this)(0,0)->get_bound(index);
//     for(int i=0; i<m;i++){
//       for(int j=0; j<n; j++){
//         ans = maximum(ans, (*this)(i,j)->get_bound(index));
//       }
//     }
//     return ans;
//   }

// };
// // compute partial derivative
//   template <unsigned int d, unsigned int m, unsigned int n, class T>
//   MVPowerseries<d,m,n,T> MVPowerseries<d,m,n,T>::partial_derivative(const unsigned int i) const{
//     std::array<unsigned int, d> index{};
//     index[i-1] = 1;
//     MVPowerseries<d,m,n,T> ans;
//     for(unsigned int i=0; i < m; i++){
//       for(unsigned int j=0; j < n; j++){
//         ans(i,j) = std::make_shared<Powerseries<d,REAL>>(std::make_shared<Cinfinity<d,REAL>>((*this)(i,j)->derive(index)), (*this)(i,j)->get_center(), (*this)(i,j)->get_radius(), (*this)(i,j)->get_M());
//       }
//     }
//     return ans;
//   }

// // redefine overloaded operators on matrix type
// template <unsigned int d, unsigned int m, unsigned int n, class T>
// MVPowerseries<d,m,n,T> operator+(const MVPowerseries<d,m,n,T>& lhs, const MVPowerseries<d,m,n,T>& rhs){
//   return add(lhs,rhs);
// }

// template <unsigned int d, unsigned int m, unsigned int n, class T>
// MVPowerseries<d,m,n,T> operator-(const MVPowerseries<d,m,n,T>& lhs, const MVPowerseries<d,m,n,T>& rhs){
//   return subtract(lhs,rhs);
// }

// template <unsigned int d, unsigned int m, unsigned int n, unsigned int k, class T>
// MVPowerseries<d,m,k,T> operator*(const MVPowerseries<d,m,n,T>& lhs, const MVPowerseries<d,n,k,T>& rhs){
//   return multiply(lhs,rhs);
// }

// template <unsigned int d, unsigned int m, unsigned int n, class T>
// MVPowerseries<d,m,n,T> operator*(const T& lhs, const MVPowerseries<d,m,n,T>& rhs){
//   MVPowerseries<d,m,n,T> ans;
//   for(int i=0; i<m;i++){
//     for(int j=0;j<n;j++){
//       ans(i,j) = lhs*rhs(i,j);
//     }
//   }
//   return ans;
// }
// template <unsigned int d, unsigned int m, unsigned int n, class T>
// MVPowerseries<d,m,n,T> operator*(const MVPowerseries<d,m,n,T>& lhs, const T& rhs){
//   return multiply(lhs,rhs);
// }

// template<unsigned int d, unsigned int m, unsigned int n,class T>
// MVPowerseries<d,m,n,T> toPowerseries(const MVFunction<d,m,n,T>& f, const std::array<T,d> center, const REAL& radius, const REAL& M){
//   MVPowerseries<d,m,n,T> ans;
//   for(int i=0; i<m;i++){
//     for(int j=0;j<n;j++){
//       ans(i,j) = std::make_shared<Powerseries<d,T>>(f(i,j), center, radius,M);
//     }
//   }
//   return ans;
// }
}


#endif
