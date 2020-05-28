/*------------------------------------------------------------
 * collection of combinatorical functions
-------------------------------------------------------------*/

#ifndef COMBINATORICS_h
#define COMBINATORICS_h
#include <vector>
#include <array>
namespace iRRAM{
  std::vector<std::vector<unsigned long>> partitions(const unsigned long n, const unsigned long k);
  std::vector<std::vector<unsigned long>> bounded_count(const std::vector<unsigned long>& bound);
  INTEGER choose(int n, int k);
  REAL inv_factorial(const unsigned int n);
  REAL inv_factorial();
  template<class... Args>
  REAL inv_factorial(const int n, Args... rest)
  {
    return inv_factorial(n)*inv_factorial(rest...);
  }
  template<unsigned int d>
  REAL inv_factorial(typename std::array<unsigned int,d>::const_iterator start){
    return inv_factorial(*start)*inv_factorial<d-1>(start+1);
  }
  template<unsigned int d>
  REAL inv_factorial(const std::array<unsigned int, d>& n)
  {
    return inv_factorial<d>(n.cbegin());
  }
  template<>
  REAL inv_factorial<1>(std::array<unsigned int,1>::const_iterator start){
    return inv_factorial(*start);
  }
  // get all partitions of size k of the number n
  std::vector<std::vector<unsigned long>> partitions(const unsigned long n, const unsigned long k){
    if(k == 1) return std::vector<std::vector<unsigned long>>{{n}};
    std::vector<std::vector<unsigned long>> ans;
    for(int i=0; i<=n; i++){
      for(std::vector<unsigned long> p : partitions(n-i, k-1)){
        p.push_back(i);
        ans.push_back(p);
      }
    }
    return ans;
  }

  INTEGER factorial(const unsigned int n){
    static std::vector<INTEGER> ans={1};
    for(int j=ans.size(); j<=n; j++){
      ans.push_back(ans.back()*j);
    }
    return ans[n];
  }

  INTEGER choose(unsigned int n, unsigned int k){
    static std::vector<std::vector<INTEGER>> mem={{1}};
    if(k > n) return 0;
    if(mem.size() > n && mem[n].size() > k)
      return mem[n][k];
    if(mem.size() <= n){
      // guarantee all needed coefficients known
      choose(n-1, k);
      mem.resize(n+1);
      mem[n] = std::vector<INTEGER>{1};
    }
    if(mem[n].size() <= k){
      choose(n, k-1);
      mem[n].resize(k+1);
    }
    if(k == 0 || n==k) mem[n][k] = 1;
    else mem[n][k] = choose(n-1,k-1)+choose(n-1,k);
    return mem[n][k];
  }

  template<unsigned int n>
  REAL choose(const std::array<unsigned int, n>& alpha, const std::array<unsigned int, n>& beta, const int size){
    if(size == 0) return 1;
    return REAL(choose(alpha[size-1],beta[size-1]))*choose<n>(alpha, beta, size-1);
  }
  template<unsigned int n>
  REAL choose(const std::array<unsigned int, n>& alpha, const std::array<unsigned int, n>& beta){
    return choose<n>(alpha, beta, n);
  }

  template<unsigned int n>
  std::vector<std::array<unsigned int, n>> bounded_count(const std::array<unsigned int,n>& bound, const int size){
    if(size == 0) return std::vector<std::array<unsigned int,n>>{std::array<unsigned int,n>()};
    std::vector<std::array<unsigned int,n>> ans; 
    auto rest=bounded_count<n>(bound, size-1);
    for(auto& v : rest){
      for(int i=0; i<=bound[size-1]; i++){
        v[size-1] = i;
        ans.push_back(v);
      }
    }
    return ans;
  }

  template<unsigned int n>
  std::vector<std::array<unsigned int, n>> bounded_count(const std::array<unsigned int,n>& bound){
    return bounded_count<n>(bound, n);
  } 
  REAL inv_factorial()
  {
    return 1;
    
  }

  template<unsigned int d>
  std::vector<std::array<unsigned int, d>> partitions(unsigned int n){
    std::vector<std::array<unsigned int,d>> ans;
    for(int i=0; i<=n;i++){
      for(auto p : partitions<d-1>(n-i)){
        std::array<unsigned int, d> partition;
        partition[0] = i;
        std::copy(p.begin(),p.end(), partition.begin()+1);
        ans.push_back(partition);
      }
    }
    return ans;
  }

  template<>
  std::vector<std::array<unsigned int, 0>> partitions<0>(unsigned int n){
    if(n > 0) return std::vector<std::array<unsigned int,0>>();
    return std::vector<std::array<unsigned int, 0>>{std::array<unsigned int, 0>()};
  }

  template<unsigned int d, class T>
T power(const std::array<T,d>& x, const std::array<unsigned int, d>& alpha){
  T ans=1;
  for(int i=0; i<d; i++){
    ans *= power(x[i], alpha[i]);
  }
  return ans;
}

REAL inv_factorial(const int n){
  using std::log;
  if ((n!=0)&&(n*log(n)-n > 2*-ACTUAL_STACK.actual_prec)){
    REAL return_value(0);
    sizetype error;
    sizetype_set(error,1,ACTUAL_STACK.actual_prec);
    return_value.seterror(error);
    return return_value;
  }
  if (n==0)
    return REAL(1);
  REAL inv_fact=inv_factorial(n-1)/REAL(n);
  return inv_fact;
}
REAL inv_factorial(const unsigned int n){
  return inv_factorial(int(n));
}
}
#endif
