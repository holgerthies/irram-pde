#ifndef FOURIER_H
#define FOURIER_H
#include "analytic.h"
#include "math.h"
namespace iRRAM{
    class Analytic_to_Fourier{
        private:
            int L;
            mutable std::vector<std::vector<std::vector<COMPLEX>>> X;
            ANALYTIC<1,REAL> f;
        public:
            COMPLEX getX(const int l, const int k, const int j) const{
                if(k < 0){
                    COMPLEX x = getX(l, -k, j);
                    return COMPLEX(real(x), -imag(x));
                }
                COMPLEX imag = COMPLEX(0,1);
                if(X.size() <= l) X.resize(l+1, std::vector<std::vector<COMPLEX>>());
                if(X[l].size() <= k) X[l].resize(k+1, std::vector<COMPLEX>());
                int J = X[l][k].size();
                if(X[l][k].size() <= j) X[l][k].resize(j+1);
                if(k == 0 && j == 0){
                    return REAL(1)/REAL(L);
                }
                if(k == 0){
                    return (REAL(1)/REAL(j+1))*power(REAL(1)/REAL(2*L), j) - power(REAL(-1)/REAL(2*L), j);
                }
                COMPLEX el = exp(-2*pi()*imag*REAL(k)*REAL(l)/REAL(L));
                COMPLEX ell = exp(-2*pi()*imag*REAL(k)*(REAL(l+1))/REAL(L));
                COMPLEX n = REAL(1)/(REAL(2)*pi()*imag*REAL(k));
                for(int nj = J; nj <= j; nj++){
                    X[l][k][nj] = n*(power(REAL(-1)/REAL(2*L), nj)*el - power(REAL(1)/REAL(2*L), nj)*ell);
                    if(nj > 0)
                        X[l][k][nj] = X[l][k][nj] + n*REAL(nj)*X[l][k][nj-1];
                }
                return X[l][k][j];
            }
        Analytic_to_Fourier(const int L, const ANALYTIC<1,REAL>& f) : L(L), f(f){}
         
        COMPLEX approx_inner(const int l, const int k, const int prec) const{
            COMPLEX ans(0,0);
            int N = -prec+1+round(log(this->f.get_M()));
            for(int j=0; j<N;j++){
                ans = ans + f.ps[l]->get_coefficient({j})*getX(l,k,j);
            }
            return ans;
        }

        COMPLEX compute_inner(const int l, const int k) const {
            REAL error = REAL(2)*(this->f.get_M());
            REAL error_factor = REAL(1)/REAL(2);
            COMPLEX ans(0,0), sum(0,0);
            sizetype trunc_error = real_to_error(error), sum_error,total_error;
            sum.geterror(sum_error);
            sizetype_add(total_error, sum_error, trunc_error);
            int j=0;
            while (sizetype_less(sum_error, trunc_error) &&
              (trunc_error.exponent >= actual_stack().actual_prec)){
                sum = sum + f.ps[l]->get_coefficient({j})*getX(l,k,j);
                error *= error_factor;
                trunc_error = real_to_error(error);
                sum.geterror(sum_error);
                sizetype curr_error;
                sizetype_add(curr_error, sum_error, trunc_error);
                if(sizetype_less(curr_error, total_error)){
                    ans = sum;
                    ans.seterror(curr_error);
                    total_error = curr_error;
                }
                j++;
            }
            return ans;
        } 
        
        COMPLEX get_coeff(const int k) const{
            COMPLEX ans(0,0);
            for(int l=0; l<L; l++){
                ans = ans + compute_inner(l,k);
            }
            return ans;
        }
    };

    class Fourier_to_Analytic{
        private:
            int L,  M;
            std::function<COMPLEX(const int)> fk;
            mutable std::vector<COMPLEX> coeffs;
        public:
        COMPLEX get_coeff(const int k) const{
            while(coeffs.size() <= k){
                coeffs.push_back(fk(k));
            }
            return coeffs[k];
        }

        Fourier_to_Analytic(const int L, const int M, const  std::function<COMPLEX(const int)>& fk) : L(L), M(M), fk(fk){}
         
        REAL get_ps_coeff(const int l, const int j) const {
            int d = 1;
            int K = 0;
            REAL error = power(REAL(2), (j+2))*REAL(M)*power(REAL(L*(d+j)), j+1)*power(REAL(1) / pi(), d)+1;
            COMPLEX ans(0,0), sum(0,0);
            sizetype trunc_error = real_to_error(error), sum_error,total_error;
            sum.geterror(sum_error);
            sizetype_add(total_error, sum_error, trunc_error);
            ans.seterror(total_error);
            COMPLEX im = COMPLEX(0,1);
            COMPLEX fac = power(2*pi(), j)*inv_factorial(j);
            if(j % 4 == 1){
                fac = fac*im;
            }
            if(j % 4 == 2){
                fac = fac*REAL(-1);
            }
            if(j % 4 == 3){
                fac = fac*(-im);
            }
            
            while (sizetype_less(sum_error, trunc_error) &&
              (trunc_error.exponent >= actual_stack().actual_prec)){
                for(int k = K; k < L*(d+j); k++){
                    if(k == 0){
                        sum = sum + fac*get_coeff(k)*exp(2*pi()*im*REAL(k)*(REAL(2*l+1)/REAL(2*L)));
                    }else{
                        sum = sum + fac*power(k,j)*get_coeff(k)*exp(-2*pi()*im*REAL(k)*(REAL(2*l+1)/REAL(2*L)));
                        COMPLEX conjk = COMPLEX(real(get_coeff(k)), -imag(get_coeff(k)));
                        sum = sum + fac*power(-k,j)*conjk*exp(2*pi()*im*REAL(k)*(REAL(2*l+1)/REAL(2*L)));
                    }
                }
                error = power(REAL(2), (j+2))*REAL(M)*power(REAL(L*(d+j)), j+1)*power(REAL(1) / pi(), d);
                trunc_error = real_to_error(error);
                sum.geterror(sum_error);
                sizetype curr_error;
                sizetype_add(curr_error, sum_error, trunc_error);
                if(sizetype_less(curr_error, total_error)){
                    ans = sum;
                    ans.seterror(curr_error);
                    total_error = curr_error;
                }
                K = L*(d+j);
                d++;
            }

            return real(ans);
        } 
        
    };
}
#endif

