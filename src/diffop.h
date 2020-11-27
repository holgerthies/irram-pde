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
  class DifferentialOperator {
  private:
    std::vector<DifferentialOperatorCoefficient<d,e,T>> coeffs;
  public:
    DifferentialOperatorCoefficient<d,e,T>& get_coefficient(const Multiindex<d>& index){
      for(auto& c : coeffs){
        if(c.get_index() == index){
          return c;
        }
      }
      coeffs.push_back(DifferentialOperatorCoefficient<d,e,T>(index, MVFunction<d,e,e,T>()));
      return coeffs[coeffs.size()-1];
    }

    DifferentialOperator(const Multiindex<d>& index, const MVFunction<d,e,e,T>& f) : coeffs{DifferentialOperatorCoefficient<d,e,T>(index,f)}{}

    DifferentialOperator(){}

    auto get_coeffs() const {return coeffs;}

    MVFunction<d,e,1,T> operator()(const MVFunction<d,e,1,T>& v) {
      MVFunction<d,e,1,T> ans;
      for(auto& coeff : coeffs){
        ans = ans+coeff(v);
      }
      return ans;
    }


  };

  template<unsigned int d, unsigned int e, class T>
  DifferentialOperator<d,e,T> operator+(const DifferentialOperator<d,e,T>& f, const DifferentialOperator<d,e,T>& g){
    DifferentialOperator<d,e,T> ans(f);
    for(auto& c : g.get_coeffs()){
      auto& coeff = ans.get_coefficient(c.get_index());
      coeff = DifferentialOperatorCoefficient<d,e,T>(coeff.get_index(),c.get_coeff()+coeff.get_coeff());
    }
    return ans;
  }
  

}
#endif
