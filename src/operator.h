#ifndef OPERATOR_H
#define OPERATOR_H
#include "cinfinity.h"

namespace iRRAM{
  template<unsigned int d, class T>
  class UnaryOperator : public Cinfinity<d,T>{
  private:
    CinfinityPtr<d,T> f;
    using coeff_fun_type = std::function<T(const CinfinityPtr<d,T>&, const Multiindex<d>&)>;
    coeff_fun_type coeff_fun;
    void update_center(const vector<T,d>& new_center) override{
      f->set_center(new_center);
      this->clear_coeffs();
    }

    T get_derivative_raw(const Multiindex<d>& index) const override{
      return coeff_fun(f, index);
    }
  public:
    UnaryOperator(const coeff_fun_type& coeff_fun, const CinfinityPtr<d,T>& f) : f(f), coeff_fun(coeff_fun) {};
    CinfinityPtr<d,T> deep_copy() override{
      return std::make_shared<UnaryOperator<d,T>>(coeff_fun, f->deep_copy());
    }
  };

  template<unsigned int d, class T>
  class BinaryOperator : public Cinfinity<d,T>{
  private:
    CinfinityPtr<d,T> lhs,rhs;
    using coeff_fun_type = std::function<T(const CinfinityPtr<d,T>&,const CinfinityPtr<d,T>&, const Multiindex<d>&)>;
    coeff_fun_type coeff_fun;
    void update_center(const vector<T,d>& new_center) override{
      lhs->set_center(new_center);
      rhs->set_center(new_center);
      this->clear_coeffs();
    }

    T get_derivative_raw(const Multiindex<d>& index) const override{
      return coeff_fun(lhs,rhs, index);
    }
  public:
    BinaryOperator(const coeff_fun_type& coeff_fun, const CinfinityPtr<d,T>& lhs, const CinfinityPtr<d,T>& rhs) : lhs(lhs), rhs(rhs), coeff_fun(coeff_fun)  {};
    CinfinityPtr<d,T> deep_copy() override{
      return std::make_shared<BinaryOperator<d,T>>(coeff_fun, lhs->deep_copy(), rhs->deep_copy());
    }
  };


  template<unsigned int d, class T>
  class ScalarOperator : public Cinfinity<d,T>{
  private:
    T lhs;
    CinfinityPtr<d,T> rhs;
    using coeff_fun_type = std::function<T(const T&,const CinfinityPtr<d,T>&, const Multiindex<d>&)>;
    coeff_fun_type coeff_fun;
    void update_center(const vector<T,d>& new_center) override{
      rhs->set_center(new_center);
      this->clear_coeffs();
    }

    T get_derivative_raw(const Multiindex<d>& index) const override{
      return coeff_fun(lhs,rhs, index);
    }
  public:
    ScalarOperator(const coeff_fun_type& coeff_fun, const T& lhs, const CinfinityPtr<d,T>& rhs) : lhs(lhs), rhs(rhs), coeff_fun(coeff_fun)  {};
    CinfinityPtr<d,T> deep_copy() override{
      return std::make_shared<ScalarOperator<d,T>>(coeff_fun, lhs, rhs->deep_copy());
    }
  };
}
#endif
