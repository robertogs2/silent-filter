/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Roberto Gutierrez
 * @Date  : 20.09.2018
 */

#ifndef SILENT_UTILITIES_HPP
#define SILENT_UTILITIES_HPP

#include <opencv2/core.hpp>    // For cv::Mat
#include <opencv2/highgui.hpp> // For cv::imread/imshow
#include "Intrinsics.hpp"
#include "MetaIntrinsics.hpp"
#include <iostream>
#include <type_traits>

namespace silent {

    /**
     *Prints a vector
     *@param inputPath
     *@param outputPath
     *@param x : vector of x positions to write
     *@param y : vector of y positions to write
     */
    inline void savePath(std::string inputPath, std::string outputPath, std::vector<size_t> &x, std::vector<size_t> &y){

    // Read the image using the OpenCV
    cv::Mat map;

    const size_t R = 200;
    const size_t G = 50;
    const size_t B = 100;

    //Loads image in Mat variable
    map = cv::imread(inputPath.c_str());

    //Both vector must have same size
    if (x.size() != y.size()){
        throw silent::Exception("Size of vectors is not equal");
    }

    cv::Vec3b newColor = cv::Vec3b(R, G, B);
    for (size_t i = 0; i < x.size(); ++i){ //Draws points in given x and y position
        map.at<cv::Vec3b>(cv::Point(y[i], x[i])) = newColor;
    }

    //Saves image
    imwrite(outputPath,map);
    }

    template<typename T>
    T bilinealInterpolation(std::vector<T> data, T xi, T yi){//Assumes points in 2d come equispaced in parallel, could be optimized 
        T eps = std::numeric_limits<T>::epsilon();
        if(std::abs(data[0] - data[1]) < eps) throw silent::Exception("Interpolation x cannot be the same");
        if(std::abs(data[2] - data[3]) < eps) throw silent::Exception("Interpolation y cannot be the same");
        T x1 = data[0]; T x2 = data[1]; T y1 = data[2]; T y2 = data[3]; 
        T f11 = data[4]; T f21 = data[5]; T f12 = data[6]; T f22 =  data[7];
        T result1 = ((xi-x2)/(x1-x2))*((yi-y2)/(y1-y2))*f11 + 
               ((xi-x1)/(x2-x1))*((yi-y2)/(y1-y2))*f21;
        T result2 =((xi-x2)/(x1-x2))*((yi-y1)/(y2-y1))*f12 +
               ((xi-x1)/(x2-x1))*((yi-y1)/(y2-y1))*f22;
               T result = result1 + result2;
        return result;
    }

    namespace fallback {

        /**
         *Prints a vector
         *@tparam T  : datatype
         *@param vec : vector to print
         */
        template<typename T>
        void printVec(std::vector<T> vec) {
            for (T val : vec) {
                std::cout << " " << val << " ";
            }
            std::cout << std::endl;
        }


        /**
         * Swaps 2 rows in a matrix, starting in an specified index of column,
         * changing the original matrix
         *
         *
         * @param[in] A a matrix
         * @param[in] row1Index index for the first row to interchange
         * @param[in] row2Index index for the second row to interchange
         * @param[in] start index of the column to start
         *
         */
        template<typename T>
        void swapRows(Matrix <T> &A, size_t row1Index, size_t row2Index, size_t start) {
            if (row1Index != row2Index) { //Checks if the swap doesn't need to happen in the same row
                for (size_t i = start;
                     i < A.cols(); ++i) { //For the entries to the right of the start column make the swap of element
                    T temp = A[row1Index][i];
                    A[row1Index][i] = A[row2Index][i];
                    A[row2Index][i] = temp;
                }
            }
        }

        /**
         * Pivots 2 rows of a matrix, so that the entry in the column of the starting row is the biggest in magnitude
         * Also it pivots a permut vector with the corresponding index
         *
         * @param[in] A a matrix to pivot
         * @param[in] columnIndex to check the entries in the rows
         * @param[in] columnStart to start changing entries in the rows
         * @param[in] rowStart to start changing entries in the rows
         * @param[in] permut vector to store the original permutation
         *
         */

        template<typename T>
        void
        pivot(Matrix <T> &A, size_t columnIndex, size_t columnStart, size_t rowStart, std::vector<size_t> &permut) {
            //Finds the maximun element in the first column to do the pivot
            T max = A[rowStart][columnIndex];
            size_t maxI = rowStart;
            for (size_t p = rowStart + 1; p < A.rows(); ++p) {
                if (std::abs(A[p][columnIndex]) > std::abs(max)) {
                    maxI = p;
                    max = A[p][columnIndex];
                }
            }
            //Swaps the row in the A matrix and in the vector
            if (maxI != rowStart) {
                std::swap(permut[rowStart], permut[maxI]);
                swapRows(A, rowStart, maxI, columnStart);
            }
        }//function pivot
    } // namespace fallback


    namespace simd {

        /**
         * Prints a matrix in a human readable way
         *
         * @param[in] A a matrix
         *
         */
        template<typename T>
        void display(Matrix <T> &A) {
            for (size_t i = 0; i < A.rows(); ++i) { //For each row in A
                for (size_t j = 0; j < A.cols(); ++j) { //For each entry in the row, print the entry
                    std::cout << A[i][j] << '\t';
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }

        /**
       * Prints a matrix in a human readable way
       *
       * @param[in] A a matrix
       *
       */
        template<typename T>
        void display(const Matrix <T> &A) {
            for (size_t i = 0; i < A.rows(); ++i) { //For each row in A
                for (size_t j = 0; j < A.cols(); ++j) { //For each entry in the row, print the entry
                    std::cout << A[i][j] << '\t';
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }

        template<typename T, typename regType>
        inline void swapRows(Matrix <T> &A, size_t row1Index, size_t row2Index) {
            regType *aptr = reinterpret_cast<regType *>(A.data());//Takes pointer to chunk of A

            const size_t tentries = A.rows() * A.dcols();//Takes amoun of entries in c
            const size_t blocks =
                    (tentries * sizeof(T) + (sizeof(regType) - 1)) / sizeof(regType); //Calculates amount of blocks
            const size_t columnBlocks = blocks / A.rows();//Calculates amount of block per row
            regType *aptr1 = aptr + row1Index * columnBlocks;
            regType *aptr2 = aptr + row2Index * columnBlocks;
            for (size_t i = 0; i < columnBlocks; ++i) {

                regType temp = *aptr1;
                *aptr1 = *aptr2;
                *aptr2 = temp;
                aptr1++;
                aptr2++;
            }


        }//function swapRows

        template<typename T, typename regType>
        inline void
        pivot(Matrix <T> &A, size_t columnIndex, size_t columnStart, size_t rowStart, std::vector<size_t> &permut) {
            //Finds the maximun element in the first column to do the pivot
            T max = A[rowStart][columnIndex];
            size_t maxI = rowStart;
            for (size_t p = rowStart + 1; p < A.rows(); ++p) {
                if (std::abs(A[p][columnIndex]) > std::abs(max)) {
                    maxI = p;
                    max = A[p][columnIndex];
                }
            }
            //Swaps the row in the A matrix and in the vector
            if (maxI != rowStart) {
                std::swap(permut[rowStart], permut[maxI]);
                swapRows<T, regType>(A, rowStart, maxI);
            }
        }//function pivot


        template<typename T, typename regType>
        regType scale(T factor, regType vector) {
            return mm_mul<T>(vector, mm_set1<T, regType>(factor));
        }

    } // namespace simd


    // The utilities implementation (uimpl) namespace
    // dispatches to the corresponding methods
#ifdef SILENT_ENABLE_SIMD
    namespace uimpl=simd;
#else
    namespace uimpl=fallback;
#endif

} // namespace silent

#endif
