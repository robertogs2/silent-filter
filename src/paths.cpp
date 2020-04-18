/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the CE3102 Numerical Analysis lecture at TEC
 *
 * @Author: 
 * @Date  : 24.02.2018
 */

#include <cstdlib>
#include <iostream>

#include <AnpiConfig.hpp>

#include <string>

#include <Matrix.hpp>
#include <Exception.hpp>

#include <DFT.hpp>

int main() {

    std::cout << "Welcome back" << std::endl;
    silent::DFT<double> D=silent::DFT<double>(4);
    std::vector<double> v = {1,2,3,4,5};
    std::vector<std::complex<double>> r;

    //r = D.transform(v);
    //r = D.transform(r, true);


    D.initOverlap({1,2});
    D.printVector(D.process(v));


    return EXIT_SUCCESS;
}
  
