/*
 * Copyright (C) 2020
 * Silent project
 * @Author: Roberto Gutierrez
 * @Date:   18.04.2020
 */


#ifndef SILENT_DFT_HPP
#define SILENT_DFT_HPP

#include <Matrix.hpp>
#include <Utilities.hpp>
#include <complex>
#include "boost/math/constants/constants.hpp"
#include <vector>

namespace silent{
    template<typename T>
    class DFT{
        public:

        size_t N;
        Matrix<std::complex<T>> W;
        Matrix<std::complex<T>> Wi;

        DFT(){
            N=0;
        }
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

        void displayW(bool inverse=false){
            Matrix<std::complex<T>> a = inverse ? Wi : W;
            printMatrix(a);
        }

        private:
        
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

    };
} 

#endif