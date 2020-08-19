#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H
#include <iRRAM.h>
#include "cinfinity.h"
#include <string>
#include <regex>
#include <sstream>
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
  class Polynomial : public Cinfinity<d,T> {
  private:
    using coeff_type = Polynomial<d-1,T>;
    std::vector<coeff_type> coeffs;
    friend class Polynomial<d+1,T>;
    // T get_coefficient_raw(const std::array<unsigned int,d>& index) const override{
    //   std::array<unsigned int, d> zero{};
    //   return get_derivative_coefficient(index.cbegin(), zero.cbegin());
    // }
    // T get_bound_raw(const std::array<unsigned int,d>& index) const override{
    //   return get_derivative_bound(index.cbegin(), this->get_radius());
    // }
    // T get_derivative_coefficient(typename std::array<unsigned int,d>::const_iterator index, typename std::array<unsigned int,d>::const_iterator derivative) const{
    //   int i = *index+*derivative;
    //   if(i > get_degree()) return 0;
    //   auto ans = FactorCache::get_factor(i,*derivative)*coeffs[i].get_derivative_coefficient(index+1, derivative+1);
    //   return ans;
    // }
    // T get_derivative_bound(typename std::array<unsigned int,d>::const_iterator index, const REAL& r) const{
    //   REAL ans = 0;
    //   for(int i=get_degree()-*index; i>=0; i--){
    //     ans = FactorCache::get_factor(*index+i,*index)*coeffs[*index+i].get_derivative_bound(index+1, r)+abs(ans)*r;
    //   }
    //   return ans;
    // }
    // REAL get_bound(const std::array<unsigned int, d>& index, const std::array<T,d>& x, const REAL& eps) const override{
    //   return get_derivative_bound(index.begin(), norm2<d>(x)+eps);
    // }
    T eval_iter(typename Multiindex<d>::const_iterator index, typename std::array<T,d>::const_iterator x) const{
      REAL ans = 0;
      for(int i=get_degree()-*index; i>=0; i--){
        ans = REAL(choose((unsigned int)(*index+i),(unsigned int)(i)))*coeffs[*index+i].eval_iter(index+1, x+1)+ans*(*x);
      }
      return ans;
    }
    T evaluate(const Multiindex<d>& index, const vector<T,d>& x) const {
      return eval_iter(index.begin(),x.M[0].begin());
    }

    T get_derivative_raw(const Multiindex<d>& index) const override{
      return evaluate(index, this->get_center());
    }

  public:


    int get_degree() const {
      return coeffs.size()-1;
    }
    Polynomial() {
    }
    Polynomial(const std::vector<coeff_type>& coeffs) : coeffs(coeffs){
    }
    CinfinityPtr<d,T> deep_copy() override{
      return std::make_shared<Polynomial<d,T>>(coeffs);
    }

    Polynomial(const std::string& s, const std::array<char,d>& variables){
      std::stringstream ss(s);
      std::string token;
      std::regex re(variables[0]+std::string("(\\^([0-9]+))?"));
      std::vector<std::string> substr;
      while (std::getline(ss, token, '+')) {
        std::smatch m; 
        if(std::regex_search(token, m, re)){
          int i=1;
          if(m[2].str().size() > 0){
            i = std::stoi(m[2]);
          }
          if(substr.size() <= i) substr.resize(i+1);
          std::string rest = m.prefix().str()+m.suffix().str();
          if(substr[i].size() == 0){
            substr[i] = rest;
          }else{
            substr[i] += '+'+rest;
          }
          if(substr[i].size() == 0) substr[i] = "1";
        }else{
          if(substr.size() == 0) substr.resize(1);
          if(substr[0].size() == 0){
            substr[0] = token;
          }else{
            substr[0] += '+'+token;
          }
        }
      }
      coeffs.resize(substr.size());
      std::array<char, d-1> vars;
      std::copy(variables.begin()+1, variables.end(), vars.begin());
      for(int i=0; i<substr.size(); i++){
        coeffs[i] = Polynomial<d-1,REAL>(substr[i], vars);
      }
    }

    T operator()(const vector<T,d>& x) const{
      Multiindex<d> zero{};
      return evaluate(zero, x);
    }

    void set_coefficient(const Multiindex<d>& index, const T& value){
      if(coeffs.size() <= index[0]){
        coeffs.resize(index[0]+1);
      }
      coeffs[index[0]].set_coefficient(remove_first(index),value);
    }
    
    void print(const std::array<char,d>& variables, const bool endline = true){
      std::array<char, d-1> vars;
      std::copy(variables.begin()+1, variables.end(), vars.begin());
      bool first = true;
      cout << "(";
      for(int i=0; i<coeffs.size(); i++){
        if(!first){
          cout << " + ";
        }
        first = false;
        coeffs[i].print(vars, false);
        if(i >= 1) cout << std::string(1,variables[0]);
        if(i > 1) cout << "^" << i;
      }
      cout << ")";
      if(endline) cout << std::endl;
    }
    template<unsigned int n,class D>
    friend Polynomial<n,D> operator*(const D& lhs, const Polynomial<n,D>& rhs);
  };

  template<class T>
  class Polynomial<0,T>  {
  private:
    T value;
    // T get_derivative_coefficient(typename std::array<unsigned int,0>::const_iterator index, typename std::array<unsigned int,0>::const_iterator derivative) const{
    //   return value;
    // }
    // T get_derivative_bound(typename std::array<unsigned int,0>::const_iterator index, const REAL& r) const{
    //   return value;
    // }
    T eval_iter(typename Multiindex<0>::const_iterator index, typename std::array<T,0>::const_iterator x) const{
      return value;
    }
    friend class Polynomial<1,T>;
  public:
    Polynomial() : value(0) {};
    Polynomial(const T& value) : value(value) {};
    Polynomial(const std::string& s, const std::array<char,0>& variables){
      std::stringstream ss(s);
      std::string token;
      while (std::getline(ss, token, '+')) {
        this->value += token;
      }
    }
    void print(const std::array<char,0>& variables, const bool endline = true){
      cout << value.as_double(2);
      if(endline) cout << std::endl;
    }
    void set_coefficient(const Multiindex<0>& index, const T& value){
      this->value = value;
    }
    template<class D>
    friend Polynomial<0,D> operator*(const D& lhs, const Polynomial<0,D>& rhs);
  };

// scalar multiplication
  template <unsigned int d, class T>
  Polynomial<d,T> operator*(const T& lhs, const Polynomial<d,T>& rhs){
    std::vector<Polynomial<0,T>> coeffs(rhs.get_degree()+1);
    for(int i=0; i<=rhs.get_degree(); i++){
      coeffs[i] = lhs*rhs.coeffs[i];
    }
    return Polynomial<d,T>(coeffs);
  }

  template <class T>
  Polynomial<0,T> operator*(const T& lhs, const Polynomial<0,T>& rhs){
    return lhs*(rhs.value);
  }
  template <unsigned int d, class T>
  Polynomial<d,T> operator*(Polynomial<d,T> const& lhs,const T& rhs){
    return rhs*lhs;
  }
  // template<unsigned int d, unsigned int m, unsigned int n,class T>
  // MVPowerseries<d,m,n,T> P2M(const std::array<std::array<std::string, n>,m>& M, const std::array<char,d>& variables){
  //   std::array<std::array<PS_ptr<d,T>,n>,m> Mptr;
  //   std::array<T, d> center{};
  //   std::array<unsigned int, d> zero{};
  //   for(int i=0; i<n; i++){
  //     for(int j=0; j<m; j++){
  //       auto poly = new Polynomial<d,T>(M[i][j], variables);
  //       std::shared_ptr<Cinfinity<d,T>> f(poly);
  //       std::shared_ptr<Powerseries<d,T>> p = std::make_shared<Powerseries<d,T>>(f,center, 1, poly->get_bound(zero,center,1)); 
  //       Mptr[i][j] = p;
  //     }
  //   }
  //   return MVPowerseries<d,m,n,T>(Mptr);
  // }
  // template<unsigned int d, unsigned int m, unsigned int n,class T>
  // MVFunction<d,m,n,T> P2MF(const std::array<std::array<std::string, n>,m>& M, const std::array<char,d>& variables){
  //   std::array<std::array<CinfinityPtr<d,T>,n>,m> Mptr;
  //   for(int i=0; i<n; i++){
  //     for(int j=0; j<m; j++){
  //       std::shared_ptr<Cinfinity<d,T>> f(new Polynomial<d,T>(M[i][j], variables));
  //       Mptr[i][j] = f;
  //     }
  //   }
  //   return MVFunction<d,m,n,T>(Mptr);
  // }
}
#endif
