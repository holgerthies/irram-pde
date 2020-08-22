#include <iRRAM.h>
#include "cinfinity.h"
#include "operator.h"
#include "funops.h"
// #include "matrix.h"
#include "powerseries.h"
#include "mfun.h"
#include "polynomial.h"
#include "diffop.h"
#include "pde.h"
using namespace iRRAM;
template<unsigned int m, unsigned int n>
void print(const Matrix<m,n,REAL>& M){
  for(int i=0; i<m; i++){
    for(int j=0;j<n;j++){
      cout << M(i,j) << " ";
    }
    cout << std::endl;
  }
}

// template<unsigned int e, unsigned int d>
// std::vector<REAL> solve(Pde_local_solution<e,d,REAL>& system, const REAL& x){
//   std::vector<REAL> sols(e);
//   auto y = system({x});
//   for(int i=0; i<e;i++){
//     sols[i] = y(i,0);
//   }
//   return sols;
// }

void compute(){
  FactorCache::init();
  std::function<REAL(const Multiindex<1>&, const vector<REAL,1>&)> f = [] (auto index, auto x) {
                                                                                           switch(index[0] % 4){
                                                                                           case 0:
                                                                                             return sin(x[0]);
                                                                                           case 1:
                                                                                             return cos(x[0]);
                                                                                           case 2:
                                                                                             return -sin(x[0]);
                                                                                           case 3:
                                                                                             return -cos(x[0]);
                                                                                           }
                                                                                           return REAL(0);
                                                                                         };
  std::function<REAL(const Multiindex<2>&, const vector<REAL,2>&)> fx = [] (auto index, auto x) {
                                                                                            if(index[1] > 0) return REAL(0);
                                                                                            switch(index[0] % 4){
                                                                                            case 0:
                                                                                              return sin(x[0]);
                                                                                            case 1:
                                                                                              return cos(x[0]);
                                                                                            case 2:
                                                                                              return -sin(x[0]);
                                                                                            case 3:
                                                                                              return -cos(x[0]);
                                                                                            }
                                                                                            return REAL(0);
                                                                                          };
  
  std::function<REAL(const Multiindex<2>&, const vector<REAL,2>&)> fy = [] (auto index, auto x) {
                                                                                            if(index[0] > 0) return REAL(0);
                                                                                            switch(index[1] % 4){
                                                                                            case 0:
                                                                                              return sin(x[1]);
                                                                                            case 1:
                                                                                              return cos(x[1]);
                                                                                            case 2:
                                                                                              return -sin(x[1]);
                                                                                            case 3:
                                                                                              return -cos(x[1]);
                                                                                            }
                                                                                            return REAL(0);
                                                                                          };
   std::shared_ptr<Cinfinity<1,REAL>> sine = std::make_shared<Cinfinity<1,REAL>>(f);
   auto cosine = derive(sine, {1});
   auto invcosine = invert(cosine);
   auto sf = REAL(2)*sine*sine+sine*REAL(3)-REAL(2)-((REAL(2)-sine/REAL(2)+REAL(2)/cosine)/(sine+cosine));
   // cout << sine->get_derivative({0}) << std::endl;
   // cout << sine->get_derivative({1}) << std::endl;
   // cout << 6*sine->get_derivative({3}) << std::endl;
   // cout << "cos" << std::endl;
   // cout << cosine->get_derivative({0}) << std::endl;
   // cout << cosine->get_derivative({1}) << std::endl;
   // cout << 2*cosine->get_derivative({2}) << std::endl;
   // cout << 6*cosine->get_derivative({3}) << std::endl;
   // cout << 24*cosine->get_derivative({4}) << std::endl;
   // cout << "s" << std::endl;
   // cout << sf->get_derivative({0}) << std::endl;
   // cout << sf->get_derivative({1}) << std::endl;
   // cout << 6*sf->get_derivative({3}) << std::endl;
   // cout << "a" << std::endl;
   // sf->set_center({REAL(0.4)});
   // cout << sf->get_derivative({0}) << std::endl;
   // cout << sf->get_derivative({1}) << std::endl;
   // cout << 6*sf->get_derivative({3}) << std::endl;
   // cout << "b" << std::endl;
   // invcosine->set_center({REAL(0.4)});
   // cout << invcosine->get_derivative({0}) << std::endl;
   // cout << invcosine->get_derivative({1}) << std::endl;
   // cout << 6*invcosine->get_derivative({3}) << std::endl;
   CinfinityPtr<2,REAL> p = std::shared_ptr<Polynomial<2,REAL>>(new Polynomial<2,REAL>("x^2+y^2",{'x','y'}));
   MVFunction<2,3,3,REAL> M1({{
       {p, p, p},
       {p, p*p*p, p*p},
       {p*p, p, p}}});
   MVFunction<2,3,3,REAL> M2({{
       {p, p-p, p},
       {p-p, p-p, p-p},
       {p-p, p-p, p}}});
   MVFunction<2,3,1,REAL> v({{{compose(cosine,p)}, {compose(sine,p)}, {p*p*p}}});
   v.set_center({1,2});
   auto vp = derive(v, {4,0});
   auto d1 = DifferentialOperator<2,3,REAL>({1,0},M1);
   auto d2 = DifferentialOperator<2,3,REAL>({0,1},M2);
   auto d = d1+d2+d2+d2+d1+d2+d2+d2;
   auto sol = Pde_solution_series<2,3,REAL>(d,v,{1,2});
   cout << sol.get_coefficient(0,4) << std::endl;
   auto dM = d(v);
   dM.set_center({0,0});
   print(dM.get_center());
   cout<<"___"<<std::endl;
   print(dM.get_derivative({0,0}));

   dM.set_center({1,2});
   cout<<"___"<<std::endl;
   print(dM.get_center());
   cout<<"___"<<std::endl;
   print(dM.get_derivative({0,0}));

   dM.set_center({2,0});
   cout<<"___"<<std::endl;
   print(dM.get_center());
   cout<<"___"<<std::endl;
   print(dM.get_derivative({0,0}));
   // auto g = compose(cosine,compose(cosine, compose(sine, compose(cosine, compose(cosine, p)))));
   // g->set_center({1,0.5});
   // auto ps = to_powerseries(compose(cosine,p), {1,2}, 1,1);
   // auto ps_trunc = ps->get_truncated_series(10);
   // ps_trunc.print({'x', 'y'});
   // cout << ps_trunc({0.2,-0.4}) << std::endl;
   // cout << ps({1.2,1.6}) << std::endl;
   // cout << ps->get_cache_size() << std::endl;
   // auto p2 = p*p;
   // p2->set_center({1,0.5});
   // cout << p->get_derivative({0,1}) << std::endl;
   // cout << p->get_derivative({1,0}) << std::endl;
   // cout << p->get_derivative({2,0}) << std::endl;
   // cout << p->get_derivative({0,2}) << std::endl;
   // cout << "comp" << std::endl;
   // cout << g->get_derivative({0,0}) << std::endl;
   // cout << g->get_derivative({1,0}) << std::endl;
   // cout << g->get_derivative({0,1}) << std::endl;
   // cout << g->get_derivative({1,1}) << std::endl;
   // cout << g->get_derivative({1,4}) << std::endl;
   // MVFunction<1,2,2,REAL> M({{{sine, compose(sine, cosine)}, {(sine+cosine)/sine, cosine}}});
   // auto M2 = derive(M,{20});
   // M2.set_center({pi()/2});
   // print(M2.get_derivative({3}));
   // cout << (*(sine+REAL(2)*sine-sine+sine-sine-sine-sine+sine*sine*(REAL(1)-(sine+sine))*(sine+sine+REAL(1))))({pi()/2}) << std::endl; 
   // auto cose = std::make_shared<Cinfinity<1,REAL>>(sine->derive({1}));
   // cout << (*cose)({{pi()}}) << std::endl; 
   // auto sc = sine*cose;
   // auto scd = sc->derive({20});
   // cout << (scd)({0.3}) << std::endl;
   // std::shared_ptr<Cinfinity<2,REAL>> sinex = std::make_shared<Cinfinity<2,REAL>>(fx, [] (auto index, auto x, auto eps) {  if(index[1] > 0) return 0; else return 1;} );
   // std::shared_ptr<Cinfinity<2,REAL>> siney = std::make_shared<Cinfinity<2,REAL>>(fy, [] (auto index, auto x, auto eps) { if(index[0] > 0) return 0; else return 1;} );
   // static int iteration_counter = 0;

   // int dimension=0, system,max_iter,prec;
   // struct rusage usage;
   // struct timeval start, end;
   // auto F1_1d = P2M<1,2,2,REAL>({
   //     {{"0", "1"},
   //      {"1", "0"},
   //     }}, {'x'});

   // auto v_1d = P2M<1,2,1,REAL>({{{"2x+0.02"}, {"1x^2+0.2"}}}, {'x'});

   // auto F1_2d = P2M<2,3,3,REAL>({
   //     {{"0", "0", "1"},
   //      {"0","0", "0"},
   //      {"1","0","0"}
   //     }}, {'x', 'y'});

   // auto F2_2d = P2M<2,3,3,REAL>({
   //     {{"0", "0", "0"},
   //      {"0","0", "1"},
   //      {"0","1","0"}
   //     }}, {'x', 'y'});

   // auto v_2d = P2M<2,3,1,REAL>({{{"2xy+0.01"}, {"1x^2+0.05"}, {"4y+0.03"}}}, {'x', 'y'});

   // auto F1_3d = P2M<3,4,4,REAL>({
   //     {{"0", "0", "0", "-1"},
   //      {"0","0", "0", "0"},
   //      {"0","0", "0", "0"},
   //      {"-1","0","0", "0"}
   //     }}, {'x', 'y', 'z'});

   // auto F2_3d = P2M<3,4,4,REAL>({
   //     {{"0", "0", "0", "0"},
   //      {"0","0", "0", "-1"},
   //      {"0","0", "0", "0"},
   //      {"0","-1","0", "0"}
   //     }}, {'x', 'y', 'z'});

   // auto F3_3d = P2M<3,4,4,REAL>({
   //     {{"0", "0", "0", "0"},
   //      {"0","0", "0", "0"},
   //      {"0","0", "0", "-1"},
   //      {"0","0","-1", "0"}
   //     }}, {'x', 'y', 'z'});

   // auto v_3d = P2M<3,4,1,REAL>({{{"2xy+0.01"}, {"1x^2+0.05"}, {"4y+0.03"}, {"3z+0.01"}}}, {'x', 'y', 'z'});
   // std::array<MVPowerseries<1,2,2,REAL>,1> farr_1d = {{F1_1d}};
   // std::array<MVPowerseries<2,3,3,REAL>,2> farr_2d = {{F1_2d, F2_2d}};

   // std::array<MVPowerseries<3,4,4,REAL>,3> farr_3d = {{F1_3d, F2_3d, F3_3d}};
   // while(dimension == 0){
   //   iRRAM::cout << "choose system dimension" << std::endl;
   //   iRRAM::cin >> dimension;
   //   if(dimension != 1 && dimension != 2 && dimension != 3){
   //     iRRAM::cout << "invalid dimension" << std::endl;
   //     dimension = 0;
   //   }
   // }
   // iRRAM::cout << "choose system" << std::endl;
   // iRRAM::cin >>  system;
  
   // REAL x;
   // iRRAM::cout << "choose x" << std::endl;
   // iRRAM::cin >>  x;

   // iRRAM::cout << "choose precision (or 0 for iteration number)" << std::endl;
   // iRRAM::cin >>  prec;
   // bool out = true;
   // if(prec > 0){
   //   iRRAM::cout << setRwidth(prec) << std::endl;
   // }
   // else{
   //   iRRAM::cout << "choose number of iterations" << std::endl;
   //   iRRAM::cin >>  max_iter;
   //   out=false;
   // }

   // iteration_counter++; 
   // if(prec <= 0 && iteration_counter == max_iter) return;
  
  
   // getrusage(RUSAGE_SELF, &usage);
   // start = usage.ru_utime;
  
   // getrusage(RUSAGE_SELF, &usage);
   // start = usage.ru_utime;
   // std::vector<REAL> ys;
   // if(dimension == 1){
   //   auto sol = Pde_local_solution<2,1,REAL>(farr_1d, v_1d, {{0}});
   //   ys = solve(sol, x);
   // }
   // if(dimension == 2){
   //   auto sol = Pde_local_solution<3,2,REAL>(farr_2d, v_2d, {{0,0}});
   //   ys = solve(sol, x);
   // }
   // if(dimension == 3){
   //   auto sol = Pde_local_solution<4,3,REAL>(farr_3d, v_3d, {{0,0,0}});
   //   ys = solve(sol, x);
   // }
   //   // for(int i = 0; i<ys.size();i++){
   //   //   cout << ys[i] << " ";
   //   // }
   //   // cout << std::endl;
   // //print(sol({x}));
   // getrusage(RUSAGE_SELF, &usage);
   // end = usage.ru_utime;
   // auto iteration_time =  end.tv_sec-start.tv_sec+double(end.tv_usec-start.tv_usec)/1000000;

   // sizetype error;
   // sizetype_exact(error);
   // for(auto y : ys){
   //   sizetype curr_error;
   //   y.geterror(curr_error);
   //   sizetype_max(error, error, curr_error);
   // }
  
   // int error_exp_normalized;
   // unsigned long mantissa = error.mantissa;
   // error_exp_normalized = error.exponent;
   // while(mantissa > 1){
   //   mantissa /= 2;
   //   error_exp_normalized++;
   // }
   // //check(sol, A3_sol(5), 50);
   // // return;
  
  
   // // std::cout << " dimension:" << dimension;
   // // std::cout << " system:" << system;
   // // std::cout << " time:" << iteration_time;
   // // std::cout << " iteration:" << iteration_counter;
   // // std::cout << " precision:" << ACTUAL_STACK.actual_prec;
   // // std::cout << " error_mantissa:"  << error.mantissa;
   // // std::cout << " error_exponent:" << error.exponent;
   // // std::cout << " normalized:" << error_exp_normalized;
   // std::cout << dimension << "," << -1*error_exp_normalized << "," << iteration_time;
  
   // std::cout << std::endl;
  
   // if(prec <= 0)
   //   iRRAM::cout << (REAL(0) == REAL(0)) << std::endl;
//print(fun1d({0.2}));
// MVPowerseries<2,2,2,REAL> fun({{{sp,cp},{cp, sp}}});
// MVPowerseries<2,2,1,REAL> v({{{sp},{cp}}});
// MVPowerseries<2,2,2,REAL> fun({{{sp,cp},{cp, sp}}});
// MVPowerseries<2,2,1,REAL> v({{{sp},{cp}}});
// // print(fun({"0", "0"}));
// // auto dv = fun.partial_derivative(1);
// // print(dv({"0.1", "0.2"}));
// // cout << std::endl;
// std::array<MVPowerseries<2,2,2,REAL>,2> farr = {{fun,fun}};
// auto sol = Pde_local_solution<2,2,REAL>(farr, v, center);
// print(sol({1}));
  
}
