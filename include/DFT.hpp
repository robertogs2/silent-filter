/*
 * Copyright (C) 2020
 * Silent project
 * @Author: Roberto Gutierrez
 * @Date:   18.04.2020
 */


#ifndef SILENT_DFT_HPP
#define SILENT_DFT_HPP

#include <Matrix.hpp>
#include <complex>
#include "boost/math/constants/constants.hpp"
#include <vector>

namespace silent{
    template<typename T>
    class DFT{
        public:
        DFT(size_t N_){
            N=N_;
            W = Matrix<std::complex<T>>(N,N);
            Wi = Matrix<std::complex<T>>(N,N);
            fillW();
            fillWi();
        }

        std::vector<std::complex<T>> 
        transform(std::vector<T> vf){
            if(vf.size() == N){
                std::vector<std::complex<T>> vc = std::vector<std::complex<T>>(vf.begin(), vf.end());
                return W*vc;
            }
            else{
                throw silent::Exception("Vector doesn't have size N");
            }
        }
        std::vector<std::complex<T>>
        transform(std::vector<std::complex<T>> vc, bool inverse=false){
            if(vc.size() == N){
                Matrix<std::complex<T>> A = inverse ? Wi : W;
                return A*vc;
            }
            else{
                throw silent::Exception("Vector doesn't have size N");
            }
        }

        void initOverlap(std::vector<T> h_filter){
            M = h_filter.size(); // Size of h(n)
            if(N <= M) throw silent::Exception("N too short for filter");
            L = N-M+1;
            h=h_filter;
            h.resize(N);
            H=transform(h);
            xm = std::vector<T>(N, T(0));
            ymc = std::vector<std::complex<T>>(N, std::complex<T>(0));
            Xm = std::vector<std::complex<T>>(N, std::complex<T>(0));
            Ym = std::vector<std::complex<T>>(N, std::complex<T>(0));
        }

        std::vector<T> process(std::vector<T> x){
            size_t L_ = 0;
            size_t input_size =x.size();
            std::vector<std::complex<T>> yc;
            while(L_ < input_size){
                if (L_ == input_size) break;
                size_t L_copy = L;
                if(L_ + L> input_size) L_copy = (input_size-L_);
                fixxm(); // Move left
                copyxm(L_copy, L_, x); // Copy right elements
                L_ += L;
                Xm = transform(xm);
                // Element wise multiplication
                std::transform(H.begin(), H.end(), 
                            Xm.begin(), Ym.begin(), 
                            std::multiplies<std::complex<T>>());
                ymc = transform(Ym, true); // Inverse transform
                yc.insert(yc.end(), ymc.begin()+(M-1), ymc.end()); // Copy elements at the end removing first M-1
            }
            std::vector<T> y;
            L_=0;
            for(std::complex<T> ele : yc){
                y.push_back(ele.real());
                if(L_>=input_size+M-2) break;
                L_++;
            }
            return y;
        }

        void displayW(bool inverse=false){
            Matrix<std::complex<T>> a = inverse ? Wi : W;
            for(size_t k=0; k<N; ++k){
                for(size_t n=0; n<N; ++n){
                    std::cout << a(k,n) << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }

        template<typename U>
        void printVector(std::vector<U> v){
            for(size_t i = 0; i < v.size(); ++i){
                std::cout << v.at(i) << ' ';
            }
            std::cout << std::endl;
        }

        private:
        size_t N;
        size_t L;
        size_t M;
        Matrix<std::complex<T>> W;
        Matrix<std::complex<T>> Wi;
        std::vector<T> h;
        std::vector<std::complex<T>> H;
        std::vector<std::complex<T>> Ym;
        std::vector<std::complex<T>> Xm;
        std::vector<T> xm;
        std::vector<std::complex<T>> ymc;

        void fillW(){
            const T pi = boost::math::constants::pi<T>();
            const T b = 2*pi/N;
            for(size_t k=0; k<N; ++k){
                for(size_t n=0;n<N;++n){
                    W(k,n)=std::polar(T(1),b*k*n);
                }
            }
        }
        void fillWi(){
            for(size_t k=0; k<N; ++k){
                for(size_t n=0;n<N;++n){
                    Wi(k,n)=std::conj(W(k,n))*T(1/T(N));
                }
            }
        }

        void fixxm(){
            std::copy(xm.begin()+L, xm.end(), xm.begin()); // Copies end to beginning
            std::fill(xm.begin()+(M-1), xm.end(),T(0)); // Sets end to 0
        }

        void copyxm(size_t L_copy, size_t L_, std::vector<T> x){
            std::copy(x.begin() + L_, x.begin() + L_ + L_copy, xm.begin()+M-1); // Check last index
        }


        
    };
} 

#endif