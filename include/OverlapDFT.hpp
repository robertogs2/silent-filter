/*
 * Copyright (C) 2020
 * Silent project
 * @Author: Roberto Gutierrez
 * @Date:   18.04.2020
 */


#ifndef SILENT_OVERLAPDFT_HPP
#define SILENT_OVERLAPDFT_HPP

#include <Matrix.hpp>
#include <DFT.hpp>

namespace silent{
    template<typename T>
    class OverlapDFT{
        public:
        DFT<T> dft;

        OverlapDFT(){
            N=0;
        }

        OverlapDFT(size_t N_){
            N = N_;
            dft = DFT<T>(N_);
        }

        OverlapDFT(DFT<T>& dft_){
            dft = dft_;
            N = dft.N;
        }

        void setFilter(std::vector<T> h_filter){
            M = h_filter.size(); // Size of h(n)
            if(N <= M) throw silent::Exception("N too short for filter");
            L = N-M+1;
            h=h_filter;
            h.resize(N);
            H=dft.transform(h);
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
                Xm = dft.transform(xm);
                // Element wise multiplication
                std::transform(H.begin(), H.end(), 
                            Xm.begin(), Ym.begin(), 
                            std::multiplies<std::complex<T>>());
                ymc = dft.transform(Ym, true); // Inverse transform
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

        private:

        size_t N;
        size_t L;
        size_t M;

        std::vector<T> h;
        std::vector<std::complex<T>> H;
        std::vector<std::complex<T>> Ym;
        std::vector<std::complex<T>> Xm;
        std::vector<T> xm;
        std::vector<std::complex<T>> ymc;

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