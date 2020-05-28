#ifndef CINFIITY_H
#define CINFIITY_H
#include <functional>
#include <iRRAM.h>
#include <combinatorics.h>
#include <matrix.h>
#include <array>
namespace iRRAM{
  template<unsigned int n, class T>
  class Cinfinity{
    using fun_type = std::function<T(const std::array<unsigned int, n>&, const std::array<T, n>&)>;
    using bound_type = std::function<REAL(const std::array<unsigned int, n>&, const std::array<T,n>&, const REAL&)>;
  protected:
    fun_type f;
    bound_type bound;
  public:
    Cinfinity() : f([] (auto index, auto x) {return 0;}), bound([] (auto index, auto x, auto eps) {return 0;}) {}
    Cinfinity(const fun_type &f, const bound_type &bound) : f(f), bound(bound) {}
    virtual ~Cinfinity() = default;
    T operator()(const std::array<T,n>&) const;
    virtual T evaluate(const std::array<unsigned int,n>& index, const std::array<T,n>& x) const{
      return f(index,x);
    }

    virtual REAL get_bound(const std::array<unsigned int, n>& index, const std::array<T,n>& x, const REAL& eps) const{
      return bound(index, x, eps);
    }

    Cinfinity<n,T> operator-();
    Cinfinity<n,T>& operator+=(const Cinfinity<n,T>& f);
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
    return evaluate(std::array<unsigned int, n>{0}, x);
  }

// derivative
template <unsigned int n, class T>
Cinfinity<n,T> Cinfinity<n,T>::derive(const std::array<unsigned int, n>& index) const{
  auto add_index = [index] (const std::array<unsigned int,n>& ind) {
                     std::array<unsigned int,n> new_index;
                     std::transform(index.begin(), index.end(), ind.begin(),new_index.begin(), std::plus<unsigned int>());
                     return new_index;
                   };
  return Cinfinity<n,T>([add_index,this] (auto ind, auto x) {return evaluate(add_index(ind), x);}, [add_index,this] (auto ind, auto x, auto eps) {return get_bound(add_index(ind), x, eps);});
}

  template <unsigned int n, class T>
  Cinfinity<n,T>& Cinfinity<n,T>::operator+=(const Cinfinity<n,T>& f){
    *this = *this + f;
    return *this;
  }

  // unary -
  template <unsigned int n, class T>
  Cinfinity<n,T> Cinfinity<n,T>::operator-(){
    return Cinfinity<n,T>([this] (auto index, auto x) {return -evaluate(index, x);}, [this] (auto ind, auto x, auto eps) {return get_bound(ind, x, eps);});
  }

// addition
    template <unsigned int n, class T>
    std::shared_ptr<Cinfinity<n,T>> operator+(std::shared_ptr<Cinfinity<n,T>> const& lhs, std::shared_ptr<Cinfinity<n,T>> const& rhs){
      return std::make_shared<Cinfinity<n,T>>([lhs,rhs] (auto index, auto x) {return lhs->evaluate(index,x)+rhs->evaluate(index,x);}, [lhs,rhs](auto index, auto x, auto eps) {return lhs->get_bound(index,x,eps)+rhs->get_bound(index,x,eps);} );
    }

// subtraction
    template <unsigned int n, class T>
    std::shared_ptr<Cinfinity<n,T>> operator-(std::shared_ptr<Cinfinity<n,T>> const& lhs, std::shared_ptr<Cinfinity<n,T>> const& rhs){
      return std::make_shared<Cinfinity<n,T>>([lhs,rhs] (auto index, auto x) {return lhs->evaluate(index,x)-rhs->evaluate(index,x);}, [lhs,rhs](auto index, auto x, auto eps) {return lhs->get_bound(index,x,eps)+rhs->get_bound(index,x,eps);} );
    }
// scalar multiplication
    template <unsigned int n, class T>
    std::shared_ptr<Cinfinity<n,T>> operator*(const T& lhs, std::shared_ptr<Cinfinity<n,T>> const& rhs){
      return std::make_shared<Cinfinity<n,T>>([lhs,rhs] (auto index, auto x) {return lhs*rhs->evaluate(index,x);}, [lhs,rhs](auto index, auto x, auto eps) {return lhs*rhs->get_bound(index,x,eps);} );
    }
  template <unsigned int n, class T>
  std::shared_ptr<Cinfinity<n,T>> operator*(std::shared_ptr<Cinfinity<n,T>> const& lhs,const T& rhs){
    return rhs*lhs;
  }
// multiplication
// a single summand of the leibniz rule
    template <unsigned int n, class T>
    T product_rule_term(const std::shared_ptr<Cinfinity<n,T>>& f, const std::shared_ptr<Cinfinity<n,T>>& g, const std::array<unsigned int,n>& alpha, const std::array<unsigned int,n>& beta, const std::array<T, n>& x){
      std::array<unsigned int,n> betap;
      std::transform(alpha.begin(), alpha.end(), beta.begin(),betap.begin(), std::minus<unsigned int>());
      return choose<n>(alpha,beta)*(f->evaluate(beta, x))*(g->evaluate(betap,x));
    }

// a bound for a single term in the leibniz rule
    template <unsigned int n, class T>
    T product_rule_term_bound(const std::shared_ptr<Cinfinity<n,T>>& f, const std::shared_ptr<Cinfinity<n,T>>& g, const std::array<unsigned int,n>& alpha, const std::array<unsigned int,n>& beta, const std::array<T, n>& x, const REAL& eps){
      std::array<unsigned int,n> betap;
      std::transform(alpha.begin(), alpha.end(), beta.begin(),betap.begin(), std::minus<unsigned int>());
      return choose<n>(alpha,beta)*(f->get_bound(beta, x, eps))*(g->get_bound(betap,x, eps));
    }

// leibniz rule

    template <unsigned int n, class T>
    T product_rule(const std::shared_ptr<Cinfinity<n,T>>& f, const std::shared_ptr<Cinfinity<n,T>>& g, const std::array<unsigned int,n>& alpha, const std::array<T, n>& x){
      T ans;
      for(auto beta : bounded_count<n>(alpha)){
        ans += product_rule_term<n>(f,g,alpha,beta,x);
      }
      return ans;
    }

    template <unsigned int n, class T>
    T product_rule_bound(const std::shared_ptr<Cinfinity<n,T>>& f, const std::shared_ptr<Cinfinity<n,T>>& g, const std::array<unsigned int,n>& alpha, const std::array<T, n>& x, const REAL& eps){
      REAL ans;
      for(auto beta : bounded_count<n>(alpha)){
        ans += product_rule_term_bound<n>(f,g,alpha,beta,x,eps);
      }
      return ans;
    }

  template <unsigned int n, class T>
  std::shared_ptr<Cinfinity<n,T>> operator*(std::shared_ptr<Cinfinity<n,T>> const& lhs, std::shared_ptr<Cinfinity<n,T>> const& rhs){
    return std::make_shared<Cinfinity<n,T>>([lhs,rhs] (auto index, auto x) {return product_rule<n>(lhs, rhs, index, x);}, [lhs,rhs](auto index, auto x, auto eps) {return product_rule_bound<n>(lhs,rhs,index,x,eps);} );
  }

// matrix valued function
    template <unsigned int d, unsigned int m, unsigned int n, class T>
      class MVFunction : public Matrix<m,n,Cinfinity<d,T>> {
    public:
      using Matrix<m,n,Cinfinity<d,T>>::Matrix;
      MVFunction(const Matrix<m,n,Cinfinity<d,T>>& M) : Matrix<m,n,Cinfinity<d,T>>(M) {};
      using Matrix<m,n,Cinfinity<d,T>>::operator();

      MVFunction<d,m,n,T> operator-();
      MVFunction<d,m,n,T> partial_derivative(const unsigned int i) const;
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
  MVFunction<d,m,n,T> MVFunction<d,m,n,T>::partial_derivative(const unsigned int i) const{
    std::array<unsigned int, d> index;
    index[i-1] = 1;
    MVFunction<d,m,n,T> ans;
    for(unsigned int i=0; i < m; i++){
      for(unsigned int j=0; j < n; j++){
        ans(i,j) = (*this)(i,j).derive(index);
      }
    }
    return ans;
  }

  // redefine overloaded operators on matrix type
  template <unsigned int d, unsigned int m, unsigned int n, class T>
  MVFunction<d,m,n,T> MVFunction<d,m,n,T>::operator-(){
    return Matrix<m,n,Cinfinity<d,T>>::operator-();
  }

  template <unsigned int d, unsigned int m, unsigned int n, class T>
  MVFunction<d,m,n,T> operator+(const MVFunction<d,m,n,T>& lhs, const MVFunction<d,m,n,T>& rhs){
    return add(lhs,rhs);
  }

  template <unsigned int d, unsigned int m, unsigned int n, class T>
  MVFunction<d,m,n,T> operator-(const MVFunction<d,m,n,T>& lhs, const MVFunction<d,m,n,T>& rhs){
    return subtract(lhs,rhs);
  }

  template <unsigned int d, unsigned int m, unsigned int n, unsigned int k, class T>
  MVFunction<d,m,k,T> operator*(const MVFunction<d,m,n,T>& lhs, const MVFunction<d,n,k,T>& rhs){
    return multiply(lhs,rhs);
  }
  template <unsigned int d, unsigned int m, unsigned int n, class T>
  MVFunction<d,m,n,T> operator*(const T& lhs, const MVFunction<d,m,n,T>& rhs){
    return multiply(lhs,rhs);
  }
  template <unsigned int d, unsigned int m, unsigned int n, class T>
  MVFunction<d,m,n,T> operator*(const MVFunction<d,m,n,T>& lhs, const T& rhs){
    return multiply(lhs,rhs);
  }
}


#endif
