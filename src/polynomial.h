#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H
#include <iRRAM.h>
#include "cinfinity.h"
namespace iRRAM{
  //forward declarations
  template<unsigned int d, class T>
  class Polynomial;

  // // upper bound of derivatives
  //   template<unsigned int d, class T>
  //   REAL upper_bound(const Polynomial<d,T>& P, const REAL& r, typename std::array<unsigned int,d>::const_iterator derivative){
  //   REAL ans = 0;
  //   for(int i=P.get_degree()-1; i>=0; i--){
  //     ans = upper_bound<d-1,T>(P.get_coeff(i), r)+abs(ans)*r;
  //   }
  //   return ans;
  // }

  // template<class T>
  // T upper_bound<0,T>(const T& P, const T& r){
  //   return abs(P);
  // }
  

  template<unsigned int d, class T>
  class Polynomial : public Powerseries<d,T> {
  private:
    using coeff_type = Polynomial<d-1,T>;
    std::vector<coeff_type> coeffs;
    friend class Polynomial<d+1,T>;
    T get_coefficient_raw(const std::array<unsigned int,d>& index) const override{
      std::array<unsigned int, d> zero{};
      return get_derivative_coefficient(index.cbegin(), zero.cbegin());
    }
    T get_bound_raw(const std::array<unsigned int,d>& index) const override{
      return get_derivative_bound(index.cbegin(), this->get_radius());
    }
    T get_derivative_coefficient(typename std::array<unsigned int,d>::const_iterator index, typename std::array<unsigned int,d>::const_iterator derivative) const{
      if(*index+*derivative > get_degree()) return 0;
      return FactorCache::get_factor(*index+*derivative,*derivative)*coeffs[*index+*derivative].get_derivative_coefficient(index+1, derivative+1);
    }
    T get_derivative_bound(typename std::array<unsigned int,d>::const_iterator index, const REAL& r) const{
      REAL ans = 0;
    for(int i=get_degree()-*index; i>=0; i--){
        ans = FactorCache::get_factor(*index+i,i)*coeffs[*index+i].get_derivative_bound(index+1, r)+abs(ans)*r;
      }
      return ans;
    }
  public:
    unsigned int get_degree() const {
      return coeffs.size()-1;
    }
    Polynomial(const std::vector<coeff_type>& coeffs, const REAL& radius) : coeffs(coeffs){
      this->radius = radius;
    }
  };
  template<class T>
  class Polynomial<0,T>  {
  private:
    T value;
    T get_derivative_coefficient(typename std::array<unsigned int,0>::const_iterator index, typename std::array<unsigned int,0>::const_iterator derivative) const{
      return value;
    }
    T get_derivative_bound(typename std::array<unsigned int,0>::const_iterator index, const REAL& r) const{
      return value;
    }
    friend class Polynomial<1,T>;
  public:
    Polynomial() : value(0) {};
    Polynomial(const T& value) : value(value) {};
  };
  template<class T>
  Polynomial<1,T> monomial(const unsigned int index, const REAL& radius){
    std::vector<Polynomial<0,T>> coeffs(index+1);
    coeffs[index] = Polynomial<0,T>(1); 
    return Polynomial<1,T>(coeffs,radius);
  }
}
#endif
