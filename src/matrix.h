#ifndef MATRIX_H
#define MATRIX_H
#include <array>
template <unsigned int m, unsigned int n, class T>
class Matrix{
private:
public:
  std::array<std::array<T, n>,m> M;
  Matrix<m,n,T>() {};
  Matrix<m,n,T>(const std::array<std::array<T,n>,m>& M) : M(M) {};
  T& operator()(const unsigned int i, const unsigned int j) {
    return M[i][j];
  } 
  
  const T& operator()(const unsigned int i, const unsigned int j) const{
    return M[i][j];
  } 

  Matrix<m,n,T> operator-();

};

// unary -
template <unsigned int m, unsigned int n, class T>
Matrix<m,n,T> Matrix<m,n,T>::operator-(){
  Matrix<m,n,T> ans;
  for(int i=0; i<m;i++){
    for(int j=0;j<n;j++){
      ans(i,j) = M(i,j);
    }
  }
  return ans;
}

// addition
  template <unsigned int m, unsigned int n, class T>
  Matrix<m,n,T> add(const Matrix<m,n,T>& lhs, const Matrix<m,n,T>& rhs){
    Matrix<m,n,T> ans;
    for(int i=0; i<m;i++){
      for(int j=0;j<n;j++){
        ans(i,j) = lhs(i,j)+rhs(i,j);
      }
    }
    return ans;
  }

// subtraction
  template <unsigned int m, unsigned int n, class T>
  Matrix<m,n,T> subtract(const Matrix<m,n,T>& lhs, const Matrix<m,n,T>& rhs){
    Matrix<m,n,T> ans;
    for(int i=0; i<m;i++){
      for(int j=0;j<n;j++){
        ans(i,j) = lhs(i,j)-rhs(i,j);
      }
    }
    return ans;
  }
// scalar multiplication
  template <unsigned int m, unsigned int n, class T>
  Matrix<m,n,T> multiply(const T& lhs, const Matrix<m,n,T>& rhs){
    Matrix<m,n,T> ans;
    for(int i=0; i<m;i++){
      for(int j=0;j<n;j++){
        ans(i,j) = lhs*rhs(i,j);
      }
    }
    return ans;
  }


// multiplication

  template <unsigned int m, unsigned int n, unsigned int k, class T>
  Matrix<m,k,T> multiply(const Matrix<m,n,T>& lhs, const Matrix<n,k,T>& rhs){
    Matrix<m,k,T> ans;
    for(unsigned int i=0; i<m;i++){
      for(unsigned int j=0;j<k;j++){
        for(unsigned int v=0; v<n; v++){
          ans(i,j) += lhs(i,v)*rhs(v,j);
        }
      }
    }
    return ans;
  }


template <unsigned int m, unsigned int n, class T>
Matrix<m,n,T> operator+(const Matrix<m,n,T>& lhs, const Matrix<m,n,T>& rhs){
  return add(lhs,rhs);
}
template <unsigned int m, unsigned int n, class T>
Matrix<m,n,T> operator-(const Matrix<m,n,T>& lhs, const Matrix<m,n,T>& rhs){
  return subtract(lhs,rhs);
}
template <unsigned int m, unsigned int n, class T>
Matrix<m,n,T> operator*(const T& lhs, const Matrix<m,n,T>& rhs){
  return multiply(lhs,rhs);
}
template <unsigned int m, unsigned int n, class T>
Matrix<m,n,T> operator*(const Matrix<m,n,T>& lhs, const T& rhs){
  return rhs*lhs;
}
template <unsigned int m, unsigned int n, unsigned int k, class T>
Matrix<m,k,T> operator*(const Matrix<m,n,T>& lhs, const Matrix<n,k,T>& rhs){
  return multiply(lhs,rhs);
}

template<class T, unsigned int n>
class vector : public Matrix<1,n,T>{
  using Matrix<1,n,T>::Matrix;
public:
  T operator[](const unsigned int& i) const{
    return (*this)(0,i);
  }
  vector(const Matrix<1,n,T>& v){
    this->M = v.M;
  }
  vector(std::initializer_list<T> l){
    int i=0;
    for(auto x : l){
      (*this)(0,i++) = x;
    }
  }
};

template<unsigned int d, class T>
T power(const vector<T,d>& x, const iRRAM::Multiindex<d>& alpha){
  T ans=1;
  for(int i=0; i<d; i++){
    ans *= power(x[i], alpha[i]);
  }
  return ans;
}

namespace iRRAM{
// euclidean norm
  template<unsigned int d>
  REAL norm2(const vector<REAL,d>& x){
    REAL ans=0;
    for(int i = 0; i<d; i++){
      ans += x[i]*x[i];
    }
    return sqrt(ans);
  }

// infinity norm
  template<unsigned int d>
  REAL inf_norm(const vector<REAL,d>& x){
    REAL ans=x[0];
    for(int i = 1; i<d; i++){
      ans = maximum(ans, x[i]);
    }
    return ans;
  }
}
#endif
