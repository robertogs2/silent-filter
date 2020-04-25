#include <Config.hpp>
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

    size_t L = 100;
    std::vector<double> vector = std::vector<double>(L);
    for(double w = 0.01; w < 3.15; w+=0.1){
        for(size_t i = 0; i < L; ++i){
            vector[i] = cos(w*i);
        }
        silent::SpectrumDFT<double> sdft = silent::SpectrumDFT<double>(L);
        double wMax = 0;
        double magMax = 0;
        sdft.findMaxW(vector, wMax, magMax);
        std::cout << w << " : " << wMax << " : " << std::abs(w-wMax) << std::endl;
    }
    
    
    return EXIT_SUCCESS;
}
  
