#ifndef MFUN_H
#define MFUN_H
#include <functional>
#include <iRRAM.h>
#include "combinatorics.h"
#include "matrix.h"
#include "powerseries.h"
#include <array>
namespace iRRAM{
// matrix valued function
  template <unsigned int d, unsigned int m, unsigned int n, class T>
  class PSMatrix : public Matrix<m,n,PS_ptr<d,T>> {
  public:
    using Matrix<m,n,PS_ptr<d,T>>::Matrix;
    PSMatrix(const Matrix<m,n,CinfinityPtr<d,T>>& M) : Matrix<m,n,CinfinityPtr<d,T>>(M) {};
    using Matrix<m,n,CinfinityPtr<d,T>>::operator();

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
