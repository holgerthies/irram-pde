# Analytic Functions and Partial Differential Equations in iRRAM

The library provides the class templates for continuous function types and PDE solving.

## Installation instructions
The library is header only, no installation is required.

The only dependency is the [iRRAM](https://github.com/fbrausse/iRRAM) framework for exact real computation.

Usage examples can be found in the `test` folder. 
To compile them the provided `Makefile` has to be adapted to your environment (e.g. the path to the irram installation has to be updated).

## Usage
The library provides PDE solvers that operate on power series and analytic function types.
The class template for power series is defined in `powerseries.h`:

```cpp
template<unsigned int d, class T>
class Powerseries { ... }
```

Which stands for a `d` dimensional power series with coefficient type `T` around some fixed center.

The simplest way to construct a power series is by directly providing the coefficient series:

```cpp
Powerseries<d,T>(const std::function<T(const Multiindex<d>&)>& coeff_fun, const vector<T,d>& center, const REAL& r, const REAL& M)
```

The `REAL` parameters `r` and `M` have to be a coefficient bound for the series.

A more practical way is to use the functional template provided in `cinfinity.h`.

```cpp
template<unsigned int d, class T>
class Cinfinity { ... }
```

`Cinfinity` is used for infinitely often differentiable functions where the sequence of derivatives around any point can be computed.
It provides the operations
```cpp
void set_center(const vector<T,d>& new_center)
T get_derivative(const Multiindex<d>& index) const{
```
`get_derivative` returns the value of the partial derivative with degree `index` at the current center which can be updated with `set_center`.
When `get_derivative` is called the derivatives are cached until the center is updated, thus enabling accessing the same derivative multiple times efficiently.

Operations on `Cinfinity` are implemented as pointer trees (using the class template `CinfinityPtr` which is a wrapper around a smart pointer to `Cinfinity`).
Currently, operators `+`, `-`, `*`, `/` are overloaded and implement the arithmetic operations.
Further function composition and a partial derivative operator are implemented.

A `Cinfinity` object can be constructed by providing the derivative function that returns partial derivatives at any point.
E.g. the `sine` function over the type `REAL` can be defined by using iRRAM's built-in `sin` and `cos` functions as follows

```cpp
Cinfinity<1,REAL> sine([] (const Multiindex<1>& m, const vector<REAL,1>& x) {
                                                 switch(k % 4){
                                                 case 0:
                                                   return sin(x[0]);
                                                 case 1:
                                                   return cos(x[9]);
                                                 case 2:
                                                   return -sin(x[0]);
                                                 case 3:default:
                                                   return cos(x[0]);
                                                 }
                                               });
```

`Cinfinity` can then be used to generate a `Powerseries` around any center by additionally providing `M` and `r`:

```cpp
template<unsigned int d, class T>
PS_ptr<d,T> to_powerseries(const CinfinityPtr<d,T>& f, const vector<T,d>& center, const REAL& r, const REAL& M)
```

The library further provides Matrix and vector classes which are used to define matrix-valued functions `MVFunction` and power series `PSMatrix` 
with point-wise operations on them.

The file `diffop.h` defines differential operators over matrix-valued function coefficients and over constant matrix coefficients.

For analytic functions on the unit hypercube `analytic.h` provides a class template

```cpp
template<unsigned int d, class T>
class ANALYTIC { ... }
```

which consists of two integer (!) parameters M and L and a d-dimensional array of powerseries each with parameters M and r=1/L that cover the hypercube.

There also is a function
```cpp
  template<unsigned int d, class T>
  ANALYTIC<d,T> to_analytic(const CinfinityPtr<d,T> f, const int L, const int M)
```
that automatically constructs an `Analytic` object from `Cinfintiy` for given parameters `L` and `M`.

### PDE solving
PDE solvers for Cauchy-Kovalevskaya type linear PDEs are contained in `pde.h`
There are two main function templates

```cpp
  template<unsigned int e, unsigned int d, class T>
  std::array<Powerseries<1,T>, e> solve_pde(const DifferentialOperator<d,e,T>& D, const MVFunction<d,e,1,REAL>& v, const vector<T,d>& x, const REAL& r, const REAL& M) 

  template<unsigned int e, unsigned int d, class T>
  std::array<Powerseries<1,T>, e> solve_pde_constant(const ConstantDifferentialOperator<d,e,T>& D, const MVFunction<d,e,1,REAL>& v, const vector<T,d>& x, const REAL& r, const REAL& M)
```

Both of them take a differential operator (which should only contain first order partial derivatives), an analytic initial value function v, a point x and parameters r and M that are a coefficient bound for 
v  and all functions that occur in the definition of the differential operator.
The second function only differs by the first one by that the Differential Operator only has constant matrix coefficients which allows a slightly faster evaluation.

For the heat equation the file `fourier.h` contains the transformation between power series and fourier series.
It is used in the implementation of the solution for the heat equation in the file `heat.cc` in the `test` folder.

## Benchmarks
For testing and evaluating the running time of our implementation, we provide some example files in the `test` subfolder.

The file `ck.cc` can be used to test the Cauchy-Kovalevskaya algorithm on several simple examples defined in `pde_examples.h`.

For each example you can choose a time t and vector x where the solution should be evaluated. 
Further, you can either choose the output precision and output the solution or choose the number of iRRAM iterations 
to evaluate the dependency of precision to running time.
In the second case no output is displayed but instead for each iteration, it outputs working precision, the resulting error bound for the solution value and the running time for the iteration.

The file `ck_double.cc` contains the same examples but does not use irram (so choosing precision and iteration number is not necessary). 
All `REAL` entities are replaced by `double` and it outputs the result in double precision.

The file `heat.cc` contains the solver for the heat equation. 
The input is similar to the one of `ck.cc`.
