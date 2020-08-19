#ifndef FUNOPS_H
#define FUNOPS_H
#include "operator.h"
namespace iRRAM{
  // addition and subtraction
  template <unsigned int d, class T>
  T add(const CinfinityPtr<d,T>& lhs, const CinfinityPtr<d,T>& rhs, const Multiindex<d>& index){
    return lhs->get_derivative(index) + rhs->get_derivative(index);
  }

  template <unsigned int d, class T>
  T subtract(const CinfinityPtr<d,T>& lhs, const CinfinityPtr<d,T>& rhs, const Multiindex<d>& index){
    return lhs->get_derivative(index) - rhs->get_derivative(index);
  }


  // scalar addition 
  template <unsigned int d, class T>
  T sadd(const T& lhs, const CinfinityPtr<d,T>& rhs, const Multiindex<d>& index){
    if(is_zero(index)){
      return rhs->get_derivative(index)+lhs;
    }
    return rhs->get_derivative(index);
  }


  // scalar multiplication
  template <unsigned int d, class T>
  T smultiply(const T& lhs, const CinfinityPtr<d,T>& rhs, const Multiindex<d>& index){
    return lhs*rhs->get_derivative(index);
  }

  // multiplication
  // a single summand of the leibniz rule
  template <unsigned int n, class T>
  T product_rule_term(const CinfinityPtr<n,T>& f, const CinfinityPtr<n,T>& g, const Multiindex<n>& alpha, const Multiindex<n>& beta){
    return (f->get_derivative(beta))*(g->get_derivative(alpha-beta));
  }


  // leibniz rule
  template <unsigned int n, class T>
  T product_rule(const CinfinityPtr<n,T>& f, const CinfinityPtr<n,T>& g, const Multiindex<n>& alpha){
    T ans;
    for(auto beta : bounded_count<n>(alpha)){
      ans += product_rule_term<n>(f,g,alpha,beta);
    }
    return ans;
  }

  template <unsigned int d, class T>
  T multiply(const CinfinityPtr<d,T>& lhs, const CinfinityPtr<d,T>& rhs, const Multiindex<d>& index){
    return product_rule(lhs, rhs, index);
  }

  // multiplicative inverse
  template<unsigned int d, class T>
  class MultiplicativeInverse : public Cinfinity<d,T>{
  private:
    CinfinityPtr<d,T> f;
    void update_center(const vector<T,d>& new_center) override{
      f->set_center(new_center);
      this->clear_coeffs();
    }

    T get_derivative_raw(const Multiindex<d>& index) const override{
      if(is_zero(index)){
        return 1/(f->get_derivative(index));
      }
      T ans;
      for(auto beta : bounded_count<d>(index)){
        if(beta != index){
          ans += this->get_derivative(beta)*f->get_derivative(index-beta);
        }
      }
      Multiindex<d> zero{};
      return (-1)*this->get_derivative(zero)*ans;
    }
  public:
    MultiplicativeInverse(const CinfinityPtr<d,T>& f) : f(f)  {};
    CinfinityPtr<d,T> deep_copy() override{
      return std::make_shared<MultiplicativeInverse<d,T>>(f->deep_copy());
    }
  };

  template <unsigned int d, class T>
  CinfinityPtr<d,T> invert(const CinfinityPtr<d,T>& f){
    return std::make_shared<MultiplicativeInverse<d,T>>(f);
  }


  // composition 

  template<unsigned int d, class T>
  class Composition : public Cinfinity<d,T>{
  private:
    CinfinityPtr<1,T> lhs;
    CinfinityPtr<d,T> rhs;
    mutable std::vector<CinfinityPtr<d,T>> rhs_pwrs;
    void update_center(const vector<T,d>& new_center) override{
      rhs->set_center(new_center);
      rhs_pwrs.clear();
      Multiindex<d> zero{};
      auto c = rhs->get_derivative(zero);
      rhs_pwrs.push_back((rhs-c));
      lhs->set_center({c});
      this->clear_coeffs();
    }
    T get_derivative_raw(const Multiindex<d>& index) const override{
      if(is_zero(index)){
        return lhs->get_derivative({0});
      }
      auto n = abs(index);
      int sz = rhs_pwrs.size();
      for(int i=sz; i<n; i++){
        rhs_pwrs.push_back(rhs_pwrs[i-1]*(rhs_pwrs[0]));
      }
      T ans = T();
      for(int i=0; i<n;i++){
        ans += lhs->get_derivative({i+1})*rhs_pwrs[i]->get_derivative(index);
      }
      return ans;
    }
  public:
    Composition(const CinfinityPtr<1,T>& lhs, const CinfinityPtr<d,T>& rhs) : lhs(lhs->deep_copy()), rhs(rhs->deep_copy()){
      this->set_center(rhs->get_center());
    }

    CinfinityPtr<d,T> deep_copy() override{
      return std::make_shared<Composition<d,T>>(lhs->deep_copy(), rhs->deep_copy());
    }

  };

  // derivative 
  template <unsigned int d, class T>
  T derivative(const Multiindex<d>& dindex, const CinfinityPtr<d,T>& f, const Multiindex<d>& index){
    return f->get_derivative(dindex+index)*get_derivative_factor(dindex,index);
  }

  template <unsigned int d, class T>
  CinfinityPtr<d,T> compose(const CinfinityPtr<1,T>& f, const CinfinityPtr<d,T>& g){
    return std::make_shared<Composition<d,T>>(f,g);
  }

  // operator overloading
  template <unsigned int d, class T>
  CinfinityPtr<d,T> operator+(const CinfinityPtr<d,T>& lhs, const CinfinityPtr<d,T>& rhs){
    return std::make_shared<BinaryOperator<d,T>>(BinaryOperator<d,T>(add<d,T>, lhs, rhs));
  }

  template <unsigned int d, class T>
  CinfinityPtr<d,T> operator-(const CinfinityPtr<d,T>& lhs, const CinfinityPtr<d,T>& rhs){
    return std::make_shared<BinaryOperator<d,T>>(BinaryOperator<d,T>(subtract<d,T>, lhs, rhs));
  }

  template <unsigned int d, class T>
  CinfinityPtr<d,T> operator*(const CinfinityPtr<d,T>& lhs, const CinfinityPtr<d,T>& rhs){
    return std::make_shared<BinaryOperator<d,T>>(BinaryOperator<d,T>(multiply<d,T>, lhs, rhs));
  }

  template <unsigned int d, class T>
  CinfinityPtr<d,T> operator*(const CinfinityPtr<d,T>& lhs, const T& rhs){
    return std::make_shared<ScalarOperator<d,T>>(ScalarOperator<d,T>(smultiply<d,T>, rhs, lhs));
  }

  template <unsigned int d, class T>
  CinfinityPtr<d,T> operator*(const T& lhs, const CinfinityPtr<d,T>& rhs){
    return std::make_shared<ScalarOperator<d,T>>(ScalarOperator<d,T>(smultiply<d,T>, lhs, rhs));
  }

  template <unsigned int d, class T>
  CinfinityPtr<d,T> operator+(const CinfinityPtr<d,T>& lhs, const T& rhs){
    return std::make_shared<ScalarOperator<d,T>>(ScalarOperator<d,T>(sadd<d,T>, rhs, lhs));
  }

  template <unsigned int d, class T>
  CinfinityPtr<d,T> operator+(const T& lhs, const CinfinityPtr<d,T>& rhs){
    return std::make_shared<ScalarOperator<d,T>>(ScalarOperator<d,T>(sadd<d,T>, lhs, rhs));
  }

  template <unsigned int d, class T>
  CinfinityPtr<d,T> operator-(const CinfinityPtr<d,T>& lhs, const T& rhs){
    return std::make_shared<ScalarOperator<d,T>>(ScalarOperator<d,T>(sadd<d,T>, (-1)*rhs, lhs));
  }

  template <unsigned int d, class T>
  CinfinityPtr<d,T> operator-(const T& lhs, const CinfinityPtr<d,T>& rhs){
    return std::make_shared<ScalarOperator<d,T>>(ScalarOperator<d,T>(sadd<d,T>, lhs, T(-1)*rhs));
  }

  template <unsigned int d, class T>
  CinfinityPtr<d,T> operator/(const CinfinityPtr<d,T>& lhs, const CinfinityPtr<d,T>& rhs){
    return lhs*invert(rhs);
  }

  template <unsigned int d, class T>
  CinfinityPtr<d,T> operator/(const T& lhs, const CinfinityPtr<d,T>& rhs){
    return lhs*invert(rhs);
  }

  template <unsigned int d, class T>
  CinfinityPtr<d,T> operator/(const CinfinityPtr<d,T>& lhs, const T& rhs){
    return (T(1)/rhs)*lhs;
  }

  template <unsigned int d, class T>
  CinfinityPtr<d,T> derive(const CinfinityPtr<d,T>& f, const Multiindex<d>& index){
    using namespace std::placeholders;
    return std::make_shared<UnaryOperator<d,T>>(UnaryOperator<d,T>(std::bind(derivative<d,T>, index, _1, _2), f));
  }
  
}

#endif
  // template <unsigned int n, class T>
  // T product_rule_bound(const CinfinityPtr<n,T>& f, const CinfinityPtr<n,T>& g, const Multiindex<n>& alpha, const vector<T,n>& x, const REAL& eps){
  //   REAL ans;
  //   for(auto beta : bounded_count<n>(alpha)){
  //     ans += product_rule_term_bound<n>(f,g,alpha,beta,x,eps);
  //   }
  //   return ans;
  // }


  // // a bound for a single term in the leibniz rule
  // template <unsigned int n, class T>
  // T product_rule_term_bound(const CinfinityPtr<n,T>& f, const CinfinityPtr<n,T>& g, const Multiindex<n>& alpha, const Multiindex<n>& beta, const vector<T, n>& x, const REAL& eps){
  //   return choose<n>(alpha,beta)*(f->get_bound(beta, x, eps))*(g->get_bound(alpha-beta,x, eps));
  // }
