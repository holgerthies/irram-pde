#ifndef MFUN_H
#define MFUN_H
#include <functional>
#include <iRRAM.h>
#include "combinatorics.h"
#include "matrix.h"
#include "cinfinity.h"
#include <array>
namespace iRRAM{
// matrix valued function
  // template<unsigned int d, class T>
  // class Cinfinity_ptr {
  // private :
  //   std::shared_ptr<Cinfinity<d,T>> f;
  // public:  
  //   Cinfinity_ptr() : f(std::make_shared<Cinfinity<d,T>>()) {};
  //   Cinfinity_ptr(const std::shared_ptr<Cinfinity<d,T>>& f) : f(f) {};
  //   Cinfinity<d,T>* operator->() const{
  //     return f.get();
  //   }
  //   T operator()(const std::array<T,d>& x) const {
  //     return (*f)(x);
  //   }

  // };
  template <unsigned int d, unsigned int m, unsigned int n, class T>
  class MVFunction : public Matrix<m,n,CinfinityPtr<d,T>> {
  public:
    using Matrix<m,n,CinfinityPtr<d,T>>::Matrix;
    MVFunction(const Matrix<m,n,CinfinityPtr<d,T>>& M) : Matrix<m,n,CinfinityPtr<d,T>>(M) {};
    using Matrix<m,n,CinfinityPtr<d,T>>::operator();

    void set_center(const vector<T,d>& new_center){
      for(int i=0; i<m; i++){
        for(int j=0; j<m; j++){
          (*this)(i,j)->set_center(new_center);
        }
      }
    }
    
    Matrix<m,n,T> get_derivative(const  Multiindex<d>& index) const{
  Matrix<m,n,T> ans;
  for(unsigned int i=0; i < m; i++){
    for(unsigned int j=0; j < n; j++){
      ans(i,j) = (*this)(i,j)->get_derivative(index);
    }
  }
  return ans;
}
MVFunction<d,m,n,T> operator-();
};


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
