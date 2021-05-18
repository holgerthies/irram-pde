#ifndef DIFFOP_H
#define DIFFOP_H
#include "cinfinity.h"
namespace iRRAM{
  template<unsigned int d, unsigned int e, class T>
  class DifferentialOperatorCoefficient{
  private:
    Multiindex<d> index;
    MVFunction<d,e,e,T> coeff;
  public:
    Multiindex<d> get_index() const{
      return index;
    }

    MVFunction<d,e,e,T> get_coeff() const{
      return coeff;
    }
    
    DifferentialOperatorCoefficient(const Multiindex<d>& index, const MVFunction<d,e,e,T>& f) : index(index), coeff(f) {};
    MVFunction<d,e,1,T> operator()(const MVFunction<d,e,1,T>& v){
      coeff.set_center(v.get_center());
      return coeff*derive(v, index);
    }
    DifferentialOperatorCoefficient(const Multiindex<d>& index) : DifferentialOperatorCoefficient(index, MVFunction<d,e,e,T>()) {};
  };

  template<unsigned int d, unsigned int e, class T>
  DifferentialOperatorCoefficient<d,e,T> operator*(const MVFunction<d,e,e,T>& g, const DifferentialOperatorCoefficient<d,e,T>& op){
    return DifferentialOperatorCoefficient<d,e,T>(op.get_index(), g*op.get_coeff());
  }

  template<unsigned int d, unsigned int e, class T>
  class DifferentialOperatorConstantCoefficient{
  private:
    Multiindex<d> index;
    Matrix<e,e,T> coeff;
  public:
    Multiindex<d> get_index() const{
      return index;
    }

    Matrix<e,e,T> get_coeff() const{
      return coeff;
    }
    
    DifferentialOperatorConstantCoefficient(const Multiindex<d>& index, const Matrix<e,e,T>& A) : index(index), coeff(A) {};
    DifferentialOperatorConstantCoefficient(const Multiindex<d>& index) : DifferentialOperatorConstantCoefficient(index, Matrix<e,e,T>()) {};

    MVFunction<d,e,1,T> operator()(const MVFunction<d,e,1,T>& v){
      return coeff*derive(v, index);
    }
  };

  template<unsigned int d, unsigned int e, class T>
  DifferentialOperatorConstantCoefficient<d,e,T> operator*(const Matrix<e,e,T>& A, const DifferentialOperatorCoefficient<d,e,T>& op){
    return DifferentialOperatorConstantCoefficient<d,e,T>(op.get_index(), A*op.get_coeff());
  }
  
  template<unsigned int d, unsigned int e, class C, class T>
  class AbstractDifferentialOperator {
  private:
    std::vector<C> coeffs;
  public:
    C& get_coefficient(const Multiindex<d>& index){
      for(auto& c : coeffs){
        if(c.get_index() == index){
          return c;
        }
      }
      coeffs.push_back(C(index));
      return coeffs[coeffs.size()-1];
    }

    template<class D>
    AbstractDifferentialOperator(const Multiindex<d>& index, const D& f) : coeffs{C(index,f)}{}

    AbstractDifferentialOperator() {}

    auto get_coeffs() const {return coeffs;}

    MVFunction<d,e,1,T> operator()(const MVFunction<d,e,1,T>& v) {
      MVFunction<d,e,1,T> ans;
      for(auto& coeff : coeffs){
        ans = ans+coeff(v);
      }
      return ans;
    }


  };

  template<unsigned int d, unsigned int e,class C, class T>
  AbstractDifferentialOperator<d,e,C,T> operator+(const AbstractDifferentialOperator<d,e,C,T>& f, const AbstractDifferentialOperator<d,e,C,T>& g){
    AbstractDifferentialOperator<d,e,C,T> ans(f);
    for(auto& c : g.get_coeffs()){
      auto& coeff = ans.get_coefficient(c.get_index());
      coeff = C(coeff.get_index(),c.get_coeff()+coeff.get_coeff());
    }
    return ans;
  }

  template<unsigned int d, unsigned int e, class T>
  using DifferentialOperator = AbstractDifferentialOperator<d,e,DifferentialOperatorCoefficient<d,e,T>,T>;

  template<unsigned int d, unsigned int e, class T>
  using ConstantDifferentialOperator = AbstractDifferentialOperator<d,e,DifferentialOperatorConstantCoefficient<d,e,T>,T>;
  

  template<unsigned int d, unsigned int e, class T>
  ConstantDifferentialOperator<d,e,T> operator*(const DifferentialOperatorConstantCoefficient<d,e,T>& c, const ConstantDifferentialOperator<d,e,T>& D){
    ConstantDifferentialOperator<d,e,T> ans;
    for(auto& di : D.get_coeffs()){
      auto& coeff = ans.get_coefficient(c.get_index()+di.get_index());
      coeff = DifferentialOperatorConstantCoefficient<d,e,T>(coeff.get_index(),c.get_coeff()*di.get_coeff());
    }
    return ans;
  }

  template<unsigned int d, unsigned int e, class T>
  ConstantDifferentialOperator<d,e,T> operator*(const ConstantDifferentialOperator<d,e,T>& lhs, const ConstantDifferentialOperator<d,e,T>& rhs){
    ConstantDifferentialOperator<d,e,T> ans;
    for(auto& di : lhs.get_coeffs()){
      ans = ans+di*rhs;
    }
    return ans;
  }
}
#endif
