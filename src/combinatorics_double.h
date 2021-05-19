/*------------------------------------------------------------
 * collection of combinatorical functions
-------------------------------------------------------------*/
#ifndef COMBINATORICS_h
#define COMBINATORICS_h
#include <vector>
#include <array>
#include <numeric>
namespace iRRAM{
  template<unsigned int d>
  class Multiindex : public std::array<int,d>{
  public:
    using std::array<int,d>::array;
    Multiindex(){
      for(int j=0; j<d; j++)
        (*this)[j] = 0;
    }
    // template<typename ...E>
    // Multiindex(E&&... e) : std::array<int,d>{{std::forward<E>(e)...}} {}
    Multiindex(std::initializer_list<int> l) {
      for(int j=0; j<d; j++)
        (*this)[j] = 0;
      int i=0;
      for(auto n : l){
        (*this)[i++] = n;
      }
    }
  };

  template<unsigned int d>
  bool is_zero(const Multiindex<d>& index){
    for(auto i : index){
      if(i != 0) return false;
    }
    return true;
  }

  template<unsigned int d>
  Multiindex<d> operator+(const Multiindex<d>& lhs, const Multiindex<d>& rhs){
    Multiindex<d> ans;
    for(int i=0; i<d; i++){
      ans[i] = lhs[i]+rhs[i];
    } 
    return ans;
  }

  template<unsigned int d>
  Multiindex<d> operator-(const Multiindex<d>& lhs, const Multiindex<d>& rhs){
    Multiindex<d> ans;
    for(int i=0; i<d; i++){
      ans[i] = lhs[i]-rhs[i];
    } 
    return ans;
  }

  template<unsigned int d>
  void print(const Multiindex<d>& idx){
    for(int i=0; i<d;i++){
      cout << idx[i] << " "; 
    }
    cout << std::endl;
  };

  template<unsigned int d>
  Multiindex<d-1> remove_first(const Multiindex<d>& index){
    Multiindex<d-1> ans;
    std::copy(index.begin()+1, index.end(), ans.begin());
    return ans;
  }

  double inv_factorial(const unsigned int n);
  double inv_factorial();

  template<unsigned int d>
  double inv_factorial_iter(typename std::array<int,d>::const_iterator start){
    return inv_factorial(*start)*inv_factorial_iter<d-1>(start+1);
  }
  template<>
  double inv_factorial_iter<1>(std::array<int,1>::const_iterator start){
    return inv_factorial(*start);
  }

  template<unsigned int d>
  double inv_factorial(const Multiindex<d>& alpha)
  {
    return inv_factorial_iter<d>(alpha.cbegin());
  }

  INTEGER factorial(const unsigned int n){
    static std::vector<INTEGER> ans={1};
    for(int j=ans.size(); j<=n; j++){
      ans.push_back(ans.back()*j);
    }
    return ans[n];
  }

  template<unsigned int d>
  double factorial(const Multiindex<d>& alpha){
    double ans =1;
    for(int i=0; i<d;i++){
      ans *= double(factorial(alpha[i]));
    }
    return ans;
  }

  template<unsigned int d>
  int abs(const Multiindex<d>& alpha){
    return std::accumulate(alpha.begin(), alpha.end(), 0);
  }

  template<unsigned int d>
  int max(const Multiindex<d>& alpha){
    return *std::max_element(alpha.begin(), alpha.end());
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
  double choose(const Multiindex<n>& alpha, const Multiindex<n>& beta, const int size){
    if(size == 0) return 1;
    return double(choose(alpha[size-1],beta[size-1]))*choose<n>(alpha, beta, size-1);
  }

  template<unsigned int n>
  double choose(const Multiindex<n>& alpha, const Multiindex<n>& beta){
    return choose<n>(alpha, beta, n);
  }

  template<unsigned int d>
  std::vector<Multiindex<d>> bounded_count(const Multiindex<d>& bound, const int size){
    if(size == 0) return std::vector<Multiindex<d>>{Multiindex<d>()};
    std::vector<Multiindex<d>> ans; 
    auto rest=bounded_count<d>(bound, size-1);
    for(auto& v : rest){
      for(int i=0; i<=bound[size-1]; i++){
        v[size-1] = i;
        ans.push_back(v);
      }
    }
    return ans;
  }

  template<unsigned int d>
  std::vector<Multiindex<d>> bounded_count(const Multiindex<d>& bound){
    return bounded_count<d>(bound, d);
  } 
  // check for deletion

  template<unsigned int d>
  double inv_factorial_iter(typename std::array<unsigned int,d>::const_iterator start){
    return inv_factorial(*start)*inv_factorial_iter<d-1>(start+1);
  }


  template<>
  double inv_factorial_iter<1>(std::array<unsigned int,1>::const_iterator start){
    return inv_factorial(*start);
  }

  template<unsigned int d>
  double inv_factorial(const std::array<unsigned int, d>& n)
  {
    return inv_factorial_iter<d>(n.cbegin());
  }
  
  std::vector<std::vector<unsigned long>> partitions(const unsigned long n, const unsigned long k);
  std::vector<std::vector<unsigned long>> bounded_count(const std::vector<unsigned long>& bound);
  INTEGER choose(int n, int k);
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



  double inv_factorial()
  {
    return 1;
    
  }

  template<unsigned int d>
  std::vector<Multiindex<d>> partitions(unsigned int n){
    std::vector<Multiindex<d>> ans;
    for(int i=0; i<=n;i++){
      for(auto p : partitions<d-1>(n-i)){
        Multiindex<d> partition;
        partition[0] = i;
        std::copy(p.begin(),p.end(), partition.begin()+1);
        ans.push_back(partition);
      }
    }
    return ans;
  }

  template<>
  std::vector<Multiindex<0>> partitions<0>(unsigned int n){
    if(n > 0) return std::vector<Multiindex<0>>();
    return std::vector<Multiindex<0>>{Multiindex<0>()};
  }

  template<unsigned int d>
  std::vector<Multiindex<d>> max_degree_indices(unsigned int n, bool with_max=true){
    std::vector<Multiindex<d>> ans;
    for(int i=0; i<=n;i++){
      bool wm=(with_max && i<n);
      for(auto p : max_degree_indices<d-1>(n, wm)){
        Multiindex<d> partition;
        partition[0] = i;
        std::copy(p.begin(),p.end(), partition.begin()+1);
        ans.push_back(partition);
      }
    }
    return ans;
  }

  template<>
  std::vector<Multiindex<0>> max_degree_indices<0>(unsigned int n, bool with_max){
    if(with_max) return std::vector<Multiindex<0>>();
    return std::vector<Multiindex<0>>{Multiindex<0>()};
  }
  


  double inv_factorial(const int n){
    using std::log;
    if (n==0)
      return double(1);
    double inv_fact=inv_factorial(n-1)/double(n);
    return inv_fact;
  }
  double inv_factorial(const unsigned int n){
    return inv_factorial(int(n));
  }

  // caches n*(n-1)*..*(n-m+1)
  class FactorCache{
  private:
    FactorCache() = default;
    ~FactorCache() = default;
    static FactorCache* instance;
    std::vector<std::vector<double>> factor_cache;
  public:
    static void init(){
      instance = new FactorCache();
    }

    double factor(int n, int m){
      while(factor_cache.size() <= n){
        factor_cache.push_back({1});
      }
      while(factor_cache[n].size() <= m){
        int j = factor_cache[n].size();
        factor_cache[n].push_back(factor_cache[n][j-1]*double(n-j+1));
      }
      return factor_cache[n][m];
    }
    static double get_factor(int n, int m){
      return instance->factor(n,m);
    }
    static auto get_size(){
      return instance->factor_cache.size();
    }
  };
  FactorCache* FactorCache::instance = 0;

  // computes (alpha+beta)!/(alpha!)
  template<unsigned int d>
  double get_derivative_factor(Multiindex<d> alpha, Multiindex<d> beta){
    double ans=1;
    for(int i=0; i<d;i++){
      ans *= FactorCache::get_factor(alpha[i]+beta[i], alpha[i]); 
    }
    return ans;
  }
}
#endif
