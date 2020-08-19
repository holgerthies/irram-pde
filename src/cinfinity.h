#ifndef CINFIITY_H
#define CINFIITY_H
#include <functional>
#include <iRRAM.h>
#include <combinatorics.h>
#include <matrix.h>
#include <array>
#include <memory>
namespace iRRAM{
  template<unsigned int d, class T> class Cinfinity;

  template<unsigned int d, class T>
  using CinfinityPtr = std::shared_ptr<Cinfinity<d,T>>;
  template<unsigned int d, class T>

  class CoefficientCache {
  private:
    std::vector<CoefficientCache<d-1,T>> cache;
  public:
    T get(typename Multiindex<d>::const_iterator start){
      return cache[*start].get(start+1);
    }

    T get(const Multiindex<d>& index){
      return get(index.cbegin());
    }

    void put(typename Multiindex<d>::const_iterator start, const T& val){
      if(cache.size() <= (*start)+1)
        cache.resize((*start)+1);
      cache[*start].put(start+1,val);
    }

    void put(const Multiindex<d>& index, const T& val){
      put(index.cbegin(), val);
    }

    bool isvalid(typename Multiindex<d>::const_iterator start){
      if(cache.size() <= *start) return false;
      return cache[*start].isvalid(start+1);
    }

    bool isvalid(const Multiindex<d>& index){
      return isvalid(index.cbegin());
    }

    void clear(){
      cache.clear();
    }

  };

  template<class T>
  class CoefficientCache<1,T> {
  private:
    std::vector<T> cache;
  public:
    T get(typename Multiindex<1>::const_iterator start){
      return cache[*start];
    }

    T get(const Multiindex<1>& index){
      return get(index.cbegin());
    }

    void put(typename Multiindex<1>::const_iterator start, const T& val){
      cache.resize(*start+1);
      cache[*start] = val;
    }

    void put(const Multiindex<1>& index, const T& val){
      put(index.cbegin(), val);
    }

    bool isvalid(typename Multiindex<1>::const_iterator start){
      return (cache.size() > *start);
    }

    bool isvalid(const Multiindex<1>& index){
      return isvalid(index.cbegin());
    }

    void clear(){
      cache.clear();
    }

  };


  template<unsigned int d, class T>
  class Cinfinity {
  public:
    using fun_type = std::function<T(const Multiindex<d>&, const vector<T,d>&)>;
  private:
    mutable CoefficientCache<d,T> coeff_cache;
    vector<T,d> center;
  protected:
    fun_type f;
    virtual T get_derivative_raw(const Multiindex<d>& index) const{
      return inv_factorial(index)*f(index,center);
    }
    virtual void update_center(const vector<T,d>& new_center){
      center = new_center;
    }
    void clear_coeffs(){
      coeff_cache.clear();
    }
  public:
    Cinfinity() : f([] (auto index, auto x) {return 0;})  {}
    Cinfinity(const fun_type &f): f(f) {}
    Cinfinity(const fun_type &f, const vector<T,d>& center): center(center), f(f) {}
    virtual ~Cinfinity() = default;
    T operator()() const;

    void set_center(const vector<T,d>& new_center){
      update_center(new_center);
      clear_coeffs();
    }

    vector<T,d> get_center() const{
      return center;
    }

    virtual CinfinityPtr<d,T> deep_copy(){
      return std::make_shared<Cinfinity<d,T>>(f,center);
    }
              

    T get_derivative(const Multiindex<d>& index) const{
      if(!coeff_cache.isvalid(index)){
        // make sure previous coefficients are already cached
        Multiindex<d> pindex(index);
        for(int i=0; i<d; i++){
          if(pindex[i] > 0){
            pindex[i]--;
            get_derivative(pindex);
            pindex[i]++;
          }
        }
        if(!coeff_cache.isvalid(index))
          coeff_cache.put(index, get_derivative_raw(index));
      }
      return coeff_cache.get(index);
    }

  };

// helper for initialization
  template<unsigned int d, class T>
  CinfinityPtr<d,T> make_cinfinity(const typename Cinfinity<d,T>::fun_type& f){
    return std::make_shared<Cinfinity<d,T>>(f);
  }

// evaluation
  template <unsigned int d, class T>
  T Cinfinity<d,T>::operator()() const{
    return get_derivative({0});
  } 


// derivative
  // template <unsigned int d, class T>
  // Cinfinity<d,T> Cinfinity<d,T>::derive(const Multiindex<d>& index) const{
  //   auto add_index = [index] (const Multiindex<d>& ind) {
  //                      return index+ind;
  //                    };
  //   return Cinfinity<d,T>([add_index,this] (auto ind, auto x) {return evaluate(add_index(ind), x);});
  // }



}

#endif
