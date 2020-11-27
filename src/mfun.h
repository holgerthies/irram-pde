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
  template <unsigned int d, unsigned int m, unsigned int n, class T>
  class MVFunction : public Matrix<m,n,CinfinityPtr<d,T>> {
  public:
    using Matrix<m,n,CinfinityPtr<d,T>>::Matrix;
    MVFunction(const Matrix<m,n,CinfinityPtr<d,T>>& M) : Matrix<m,n,CinfinityPtr<d,T>>(M) {};
    MVFunction() {
      for(int i=0; i<m;i++){
        for(int j=0; j<n; j++){
          (*this)(i,j) = std::make_shared<Cinfinity<d,T>>();
        }
      }
    };
    using Matrix<m,n,CinfinityPtr<d,T>>::operator();

    void set_center(const vector<T,d>& new_center){
      for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
          (*this)(i,j)->set_center(new_center);
        }
      }
    }

    vector<T,d> get_center() const{
      return (*this)(0,0)->get_center();
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
    MVFunction<d,m,n,T>& operator+=(const MVFunction<d,m,n,T>& f){
      *this = *this + f;
      return *this;
    }
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
  return multiply(MVFunction<d,m,k,T>(), lhs,rhs);
}

template <unsigned int d, unsigned int m, unsigned int n, class T>
MVFunction<d,m,n,T> operator*(const T& lhs, const MVFunction<d,m,n,T>& rhs){
  MVFunction<d,m,n,T> ans;
    for(int i=0; i<m;i++){
      for(int j=0;j<n;j++){
        ans(i,j) = lhs*rhs(i,j);
      }
    }
    return ans;
}

template <unsigned int d, unsigned int m, unsigned int n, class T>
MVFunction<d,m,n,T> operator*(const MVFunction<d,m,n,T>& lhs, const T& rhs){
  return multiply(lhs,rhs);
}

template <unsigned int d, unsigned int e, class T>
MVFunction<d,e,e,T> identity(){
  MVFunction<d,e,e,T> ans;
  for(int i=0; i<e; i++){
    ans(i,i) = std::make_shared<Cinfinity<d,T>>(REAL(1));
  }
  return ans;
}

  template <unsigned int d, unsigned int m, unsigned int n, class T>
  MVFunction<d,m,n,T> derive(const MVFunction<d,m,n,T>& M, const  Multiindex<d>& index){
    MVFunction<d,m,n,T> ans;
    for(unsigned int i=0; i < m; i++){
      for(unsigned int j=0; j < n; j++){
        ans(i,j) = derive(M(i,j), index);
      }
    }
    return ans;
  }
}


#endif
