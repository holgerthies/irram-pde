#ifndef CINFIITY_H
#define CINFIITY_H
#include <functional>
#include <iRRAM.h>
#include <combinatorics.h>
#include <matrix.h>
#include <array>
#include <memory>
namespace iRRAM{

template<unsigned int d, class T> class Cinfinity;

template<unsigned int d, class T>
using CinfinityPtr = std::shared_ptr<Cinfinity<d,T>>;

template<unsigned int d, class T>
class Cinfinity : public std::enable_shared_from_this<Cinfinity<d,T>>{
  using fun_type = std::function<T(const Multiindex<d>&, const vector<T,d>&)>;
  using bound_type = std::function<REAL(const Multiindex<d>&, const vector<T,d>&, const REAL&)>;
protected:
  fun_type f;
  bound_type bound;
public:
  Cinfinity() : f([] (auto index, auto x) {return 0;}), bound([] (auto index, auto x, auto eps) {return 0;}) {}
  Cinfinity(const fun_type &f, const bound_type &bound) : f(f), bound(bound) {}
  virtual ~Cinfinity() = default;
  T operator()(const vector<T,d>&) const;
  virtual T evaluate(const Multiindex<d>& index, const vector<T,d>& x) const{
    return f(index,x);
  }

  virtual REAL get_bound(const Multiindex<d>& index, const vector<T,d>& x, const REAL& eps) const{
    return bound(index, x, eps);
  }
  virtual Cinfinity<d,T> add(const CinfinityPtr<d,T>& f) const;
  virtual Cinfinity<d,T> subtract(const CinfinityPtr<d,T>& f) const;
  virtual Cinfinity<d,T> multiply(const CinfinityPtr<d,T>& f) const;
  virtual Cinfinity<d,T> multiply(const T& m) const;
  Cinfinity<d,T> derive(const Multiindex<d>& index) const;
};

// evaluation
template <unsigned int d, class T>
T Cinfinity<d,T>::operator()(const vector<T,d>& x) const{
  return evaluate({0}, x);
} 

template <unsigned int n, class T>
Cinfinity<n,T> Cinfinity<n,T>::add(const CinfinityPtr<n,T>& f) const{
  auto self = this->shared_from_this();
  return Cinfinity<n,T>([self,f] (auto index, auto x) {return self->evaluate(index,x)+f->evaluate(index,x);}, [self,f](auto index, auto x, auto eps) {return self->get_bound(index,x,eps)+f->get_bound(index,x,eps);} );
}

template <unsigned int n, class T>
Cinfinity<n,T> Cinfinity<n,T>::subtract(const CinfinityPtr<n,T>& f) const{
  auto self = this->shared_from_this();
  return Cinfinity<n,T>([self,f] (auto index, auto x) {return self->evaluate(index,x)-f->evaluate(index,x);}, [self,f](auto index, auto x, auto eps) {return self->get_bound(index,x,eps)+f->get_bound(index,x,eps);} );
}

// a single summand of the leibniz rule
template <unsigned int n, class T>
T product_rule_term(const Cinfinity<n,T>& f, const Cinfinity<n,T>& g, const Multiindex<n>& alpha, const Multiindex<n>& beta, const vector<T, n>& x){
  return choose<n>(alpha,beta)*(f.evaluate(beta, x))*(g.evaluate(alpha-beta,x));
}

// a bound for a single term in the leibniz rule
template <unsigned int n, class T>
T product_rule_term_bound(const Cinfinity<n,T>& f, const Cinfinity<n,T>& g, const Multiindex<n>& alpha, const Multiindex<n>& beta, const vector<T, n>& x, const REAL& eps){
  return choose<n>(alpha,beta)*(f.get_bound(beta, x, eps))*(g.get_bound(alpha-beta,x, eps));
}

// leibniz rule

template <unsigned int n, class T>
T product_rule(const Cinfinity<n,T>& f, const Cinfinity<n,T>& g, const Multiindex<n>& alpha, const vector<T,n>& x){
  T ans;
  for(auto beta : bounded_count<n>(alpha)){
    ans += product_rule_term<n>(f,g,alpha,beta,x);
  }
  return ans;
}

template <unsigned int n, class T>
T product_rule_bound(const Cinfinity<n,T>& f, const Cinfinity<n,T>& g, const Multiindex<n>& alpha, const vector<T,n>& x, const REAL& eps){
  REAL ans;
  for(auto beta : bounded_count<n>(alpha)){
    ans += product_rule_term_bound<n>(f,g,alpha,beta,x,eps);
  }
  return ans;
}

template <unsigned int n, class T>
Cinfinity<n,T> Cinfinity<n,T>::multiply(const CinfinityPtr<n,T>& f) const{
  auto self = this->shared_from_this();
  return Cinfinity<n,T>([self,f] (const Multiindex<n>& index, auto x) {return product_rule<n>(*self, *f, index, x);}, [self,f](const Multiindex<n>& index, auto x, auto eps) {return product_rule_bound<n>(*self,*f,index,x,eps);} );
}

template <unsigned int n, class T>
Cinfinity<n,T> Cinfinity<n,T>::multiply(const T& m) const{
  auto self = this->shared_from_this();
  return Cinfinity<n,T>([self,m] (auto index, auto x) {return m*self->evaluate(index,x);}, [self,m](auto index, auto x, auto eps) {return m*self->get_bound(index,x,eps);} );
}

// derivative
template <unsigned int d, class T>
Cinfinity<d,T> Cinfinity<d,T>::derive(const Multiindex<d>& index) const{
  auto add_index = [index] (const Multiindex<d>& ind) {
                     return index+ind;
                   };
  return Cinfinity<d,T>([add_index,this] (auto ind, auto x) {return evaluate(add_index(ind), x);}, [add_index,this] (auto ind, auto x, auto eps) {return get_bound(add_index(ind), x, eps);});
}

// addition
template <unsigned int n, class T>
CinfinityPtr<n,T> operator+(const CinfinityPtr<n,T>& lhs, const CinfinityPtr<n,T>& rhs){
  return std::make_shared<Cinfinity<n,T>>(lhs->add(rhs));
}

// subtraction
  template <unsigned int n, class T>
  std::shared_ptr<Cinfinity<n,T>> operator-(std::shared_ptr<Cinfinity<n,T>> const& lhs, std::shared_ptr<Cinfinity<n,T>> const& rhs){
    return std::make_shared<Cinfinity<n,T>>(lhs->subtract(rhs));
  }
// scalar multiplication
template <unsigned int n, class T>
std::shared_ptr<Cinfinity<n,T>> operator*(const T& lhs, const CinfinityPtr<n,T>& rhs){
  return std::make_shared<Cinfinity<n,T>>(rhs->multiply(lhs));
}
template <unsigned int n, class T>
std::shared_ptr<Cinfinity<n,T>> operator*(std::shared_ptr<Cinfinity<n,T>> const& lhs,const T& rhs){
  return rhs*lhs;
}
// // multiplication

  template <unsigned int n, class T>
  std::shared_ptr<Cinfinity<n,T>> operator*(std::shared_ptr<Cinfinity<n,T>> const& lhs, std::shared_ptr<Cinfinity<n,T>> const& rhs){
    return std::make_shared<Cinfinity<n,T>>(lhs->multiply(rhs));
  }

//   template<unsigned int n,class T>
//   std::shared_ptr<Cinfinity<n,T>> from1d(const std::shared_ptr<Cinfinity<1,T>>& f, const unsigned int i){
//     return std::make_shared<Cinfinity<n,T>>(
//       [f,i] (const std::array<unsigned int, n>& index, const std::array<T, n>& x) {
//         for(int j=0; j<n;j++){
//           if(i != j && index[j] > 0) return T();
//         }
//         return f->evaluate({index[i]}, {x[i]});
//       }, 
//       [f,i] (const std::array<unsigned int, n>& index, const std::array<T, n>& x, const REAL& eps) {
//         for(int j=0; j<n;j++){
//           if(i != j && index[j] > 0) return T();
//         }
//         return f->get_bound({index[i]}, {x[i]}, eps);
//       } 
//       ); 
//   }
}

#endif
