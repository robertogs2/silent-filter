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

    template<typename U>
    void printVector(const std::vector<U> v){
        for(size_t i = 0; i < v.size(); ++i){
            std::cout << v.at(i) << ' ';
        }
        std::cout << std::endl;
    }

    template<typename U>
    void printMatrix(const Matrix<U> M){
        for(size_t k=0; k<M.rows(); ++k){
            for(size_t n=0; n<M.cols(); ++n){
                std::cout << M(k,n) << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    namespace fallback {
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
