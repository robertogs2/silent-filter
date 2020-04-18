/*
 * Copyright (C) 2020
 * Silent project
 * @Author: Roberto Gutierrez
 * @Date:   18.04.2020
 */


#ifndef SILENT_SPRECTRUMDFT_HPP
#define SILENT_SPRECTRUMDFT_HPP

#include <Matrix.hpp>
#include <DFT.hpp>

namespace silent{
    template<typename T>
    class SpectrumDFT{
        public:
        DFT<T> dft;

        SpectrumDFT(){
            N=0;
        }

        SpectrumDFT(size_t N_){
            N = N_;
            dft = DFT<T>(N_);
        }

        SpectrumDFT(DFT<T>& dft_){
            dft = dft_;
            N = dft.N;
        }

        void findMaxK(std::vector<T> x, T& wMax, T& magMax){
            if (x.size() > N) throw silent::Exception("N too short for spectrum");
            std::vector<std::complex<T>> Xk = dft.transform(x);
            const T pi = boost::math::constants::pi<T>();
            size_t k = 0;
            size_t kMax = 0;
            T mag= 0;
            magMax = 0;
            for(std::complex<T> ele : Xk){
                mag = std::abs(ele);
                if(mag > magMax){
                    magMax = mag;
                    kMax = k;
                }
                ++k; 
            }
            wMax = kMax*2*pi/T(N);

        }

        private:

        size_t N;

    };
}

#endif