/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the CE3102 Numerical Analysis lecture at TEC
 *
 * @Author: 
 * @Date  : 24.02.2018
 */

#include <AnpiConfig.hpp>
#include <Matrix.hpp>
#include <Exception.hpp>
#include <DFT.hpp>
#include <OverlapDFT.hpp>
#include <SpectrumDFT.hpp>
#include <Utilities.hpp>

#include <cstdlib>
#include <iostream>
#include <string>


int main() {

    std::cout << "Welcome back" << std::endl;
    silent::DFT<double> D=silent::DFT<double>(4);
    std::vector<double> v = {1,2,3,4,5};
    std::vector<std::complex<double>> r;
    silent::OverlapDFT<double> odft = silent::OverlapDFT<double>(D);

    //r = D.transform(v);
    //r = D.transform(r, true);


    // odft.setFilter({1,2});
    // silent::printVector(odft.process(v));

    std::vector<double> vector = std::vector<double>(1000);
    for(size_t i = 0; i < 1000; ++i){
        vector[i] = cos(0.2*i);
    }
    silent::SpectrumDFT<double> sdft = silent::SpectrumDFT<double>(1000);
    double wMax = 0;
    double magMax = 0;
    sdft.findMaxK(vector, wMax, magMax);
    std::cout << wMax << std::endl;
    return EXIT_SUCCESS;
}
  
