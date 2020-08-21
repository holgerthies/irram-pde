#ifndef DIFFOP_H
#define DIFFOP_H
#include <iRRAM.h>
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
  };

  template<unsigned int d, unsigned int e, class T>
  DifferentialOperatorCoefficient<d,e,T> operator*(const MVFunction<d,e,e,T>& g, const DifferentialOperatorCoefficient<d,e,T>& op){
    return DifferentialOperatorCoefficient<d,e,T>(op.get_index(), g*op.get_coeff());
  }

  
  template<unsigned int d, unsigned int e, class T>
  DifferentialOperatorCoefficient<d,e,T> diff(const int var, const int power){
    Multiindex<d> index{};
    index[var-1] = power;
    return DifferentialOperatorCoefficient<d,e,T>(index, identity<d,e,T>());
  }
   

  template<unsigned int d, unsigned int e, class T>
  class DifferentialOperator {
  private:
    std::vector<DifferentialOperatorCoefficient<d,e,T>> coeffs;
  public:
    MVFunction<d,e,1,T> operator()(const MVFunction<d,e,1,T>& v) const{
      MVFunction<d,e,1,T> ans;
      for(auto& coeff : coeffs){
        ans += coeff(v);
      }
      return ans;
    }

  };

}
#endif
