/*
 * Copyright (C) 2017 
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Pablo Alvarado
 * @Date:   28.12.2017
 */

#ifndef SILENT_MATRIX_ARITHMETIC_HPP
#define SILENT_MATRIX_ARITHMETIC_HPP

#include "Intrinsics.hpp"
#include "MetaIntrinsics.hpp"
#include "Utilities.hpp"
#include <iostream>
#include <type_traits>

namespace silent{
  namespace fallback {
    /*
     * Sum
     */

    // Fallback implementation
    
    // In-copy implementation c=a+b
    template<typename T,class Alloc>
    inline void add(const Matrix<T,Alloc>& a,
                    const Matrix<T,Alloc>& b,
                    Matrix<T,Alloc>& c) {

      assert( (a.rows() == b.rows()) &&
              (a.cols() == b.cols()) );

      const size_t tentries = a.rows()*a.dcols();
      c.allocate(a.rows(),a.cols());
      
      T* here        = c.data();
      T *const end   = here + tentries;
      const T* aptr = a.data();
      const T* bptr = b.data();

      for (;here!=end;) {
        *here++ = *aptr++ + *bptr++;
      }
    }

    // In-place implementation a = a+b
    template<typename T,class Alloc>
    inline void add(Matrix<T,Alloc>& a,
                    const Matrix<T,Alloc>& b) {

      assert( (a.rows() == b.rows()) &&
              (a.cols() == b.cols()) );

      const size_t tentries = a.rows()*a.dcols();
      
      T* here        = a.data();
      T *const end   = here + tentries;
      
      const T* bptr = b.data();

      for (;here!=end;) {
        *here++ += *bptr++;
      }
    }


    /*
     * Subtraction
     */

    // Fall back implementations

    // In-copy implementation c=a-b
    template<typename T,class Alloc>
    inline void subtract(const Matrix<T,Alloc>& a,
                         const Matrix<T,Alloc>& b,
                         Matrix<T,Alloc>& c) {

      assert( (a.rows() == b.rows()) &&
              (a.cols() == b.cols()) );

      const size_t tentries = a.rows()*a.dcols();
      c.allocate(a.rows(),a.cols());
      
      T* here        = c.data();
      T *const end   = here + tentries;
      const T* aptr = a.data();
      const T* bptr = b.data();

      for (;here!=end;) {
        *here++ = *aptr++ - *bptr++;
      }
    }

    // In-place implementation a = a-b
    template<typename T,class Alloc>
    inline void subtract(Matrix<T,Alloc>& a,
                         const Matrix<T,Alloc>& b) {

      assert( (a.rows() == b.rows()) &&
              (a.cols() == b.cols()) );
      
      const size_t tentries = a.rows()*a.dcols();
      
      T* here        = a.data();
      T *const end   = here + tentries;
      
      const T* bptr = b.data();

      for (;here!=end;) {
        *here++ -= *bptr++;
      }
    }

    /*
     * Matrix Multiplication
     */

    // In-place implementation a = a*b
    template<typename T,class Alloc>
    inline void multiply(Matrix<T,Alloc>& a,
                         const Matrix<T,Alloc>& b) {
      assert(a.cols() == b.rows());
      Matrix<T,Alloc> c;
      c.allocate(a.rows(),b.cols());
      for(size_t i = 0; i < a.rows(); ++i){//For each column in a
        for(size_t j = 0; j < b.cols(); ++j){ //For each column in b
          c[i][j] = T(0);
          for(size_t k = 0; k < a.cols(); ++k){
            c[i][j] += a[i][k] * b[k][j]; 
          }
        }    
      }
      a = c;
    }//function multiply

    // In-copy implementation c = a*b
    template<typename T,class Alloc>
    inline void multiply(const Matrix<T,Alloc>& a,
                         const Matrix<T,Alloc>& b,
                         Matrix<T,Alloc>& c) {
      assert(a.cols() == b.rows());
      c.allocate(a.rows(),b.cols());
      for(size_t i = 0; i < a.rows(); ++i){//For each column in a
        for(size_t j = 0; j < b.cols(); ++j){ //For each column in b
          c[i][j] = T(0);
          for(size_t k = 0; k < a.cols(); ++k){
            c[i][j] += a[i][k] * b[k][j]; 
          }
        }    
      }
    }//function multiply


    /*
     * Scalar Multiplication
     */

    // In-place implementation a = a*b
    template<typename T,class Alloc>
    inline void multiply(Matrix<T,Alloc>& a,
                         const T& b) {
      
      const size_t tentries = a.rows()*a.dcols();
      
      T* here        = a.data();
      T *const end   = here + tentries;
      

      for (;here!=end;) {
        *here++ *= b;
      }
    }//function multiply

    // In-copy implementation c = a*b
    template<typename T,class Alloc>
    inline void multiply(const Matrix<T,Alloc>& a,
                         const T& b,
                         Matrix<T,Alloc>& c) {

      const size_t tentries = a.rows()*a.dcols();
      c.allocate(a.rows(),a.cols());
      
      T* here        = c.data();
      T *const end   = here + tentries;
      const T* aptr = a.data();

      for (;here!=end;) {
        *here++ = *aptr++ * b;
      }
    }//function multiply
  } // namespace fallback

  namespace simd
  {
    
    /*
     * Addition
    */
    
    // On-copy implementation c=a+b
    template<typename T,class Alloc,typename regType>
    inline void addSIMD(const Matrix<T,Alloc>& a, 
                        const Matrix<T,Alloc>& b,
                        Matrix<T,Alloc>& c) {

      // This method is instantiated with unaligned allocators.  We
      // allow the instantiation although externally this is never
      // called unaligned
      static_assert(!extract_alignment<Alloc>::aligned ||
		    (extract_alignment<Alloc>::value >= sizeof(regType)),
		    "Insufficient alignment for the registers used");
      
      const size_t tentries = a.rows()*a.dcols();
      c.allocate(a.rows(),a.cols());

      regType* here        = reinterpret_cast<regType*>(c.data());
      const size_t  blocks = ( tentries*sizeof(T) + (sizeof(regType)-1) )/
        sizeof(regType);
      regType *const end   = here + blocks;
      const regType* aptr  = reinterpret_cast<const regType*>(a.data());
      const regType* bptr  = reinterpret_cast<const regType*>(b.data());
      
      for (;here!=end;) {
        *here++ = mm_add<T>(*aptr++,*bptr++);
      }
      
    }
       
    // On-copy implementation c=a+b for SIMD-capable types
    template<typename T,
	     class Alloc,
	     typename std::enable_if<is_simd_type<T>::value,int>::type=0>
    inline void add(const Matrix<T,Alloc>& a,
                    const Matrix<T,Alloc>& b,
                    Matrix<T,Alloc>& c) {

      assert( (a.rows() == b.rows()) &&
              (a.cols() == b.cols()) );


      if (is_aligned_alloc<Alloc>::value) {        
#ifdef __AVX512F__
        addSIMD<T,Alloc,typename avx512_traits<T>::reg_type>(a,b,c);
#elif  __AVX__
        addSIMD<T,Alloc,typename avx_traits<T>::reg_type>(a,b,c);
#elif  __SSE2__
        addSIMD<T,Alloc,typename sse2_traits<T>::reg_type>(a,b,c);
#else
        ::silent::fallback::add(a,b,c);
#endif
      } else { // allocator seems to be unaligned
        ::silent::fallback::add(a,b,c);
      }
    }

    // Non-SIMD types such as complex
    template<typename T,
             class Alloc,
             typename std::enable_if<!is_simd_type<T>::value,int>::type = 0>
    inline void add(const Matrix<T,Alloc>& a,
                    const Matrix<T,Alloc>& b,
                    Matrix<T,Alloc>& c) {
      
      ::silent::fallback::add(a,b,c);
    }

    // In-place implementation a = a+b
    template<typename T,class Alloc>
    inline void add(Matrix<T,Alloc>& a,
                    const Matrix<T,Alloc>& b) {

      add(a,b,a);
    }

    /*
     * Substraction
    */

    // On-copy implementation c=a-b
    template<typename T,class Alloc,typename regType>
    inline void subSIMD(const Matrix<T,Alloc>& a, 
                        const Matrix<T,Alloc>& b,
                        Matrix<T,Alloc>& c) {

      // This method is instantiated with unaligned allocators.  We
      // allow the instantiation although externally this is never
      // called unaligned
      static_assert(!extract_alignment<Alloc>::aligned ||
        (extract_alignment<Alloc>::value >= sizeof(regType)),
        "Insufficient alignment for the registers used");
      const size_t tentries = a.rows()*a.dcols();
      c.allocate(a.rows(),a.cols());

      regType* here        = reinterpret_cast<regType*>(c.data());
      const size_t  blocks = ( tentries*sizeof(T) + (sizeof(regType)-1) )/
        sizeof(regType);
      regType *const end   = here + blocks;
      const regType* aptr  = reinterpret_cast<const regType*>(a.data());
      const regType* bptr  = reinterpret_cast<const regType*>(b.data());
      
      for (;here!=end;) {
        *here++ = mm_sub<T>(*aptr++,*bptr++);
      }
      
    }

    // In-copy implementation c=a-b
    template<typename T,class Alloc,
       typename std::enable_if<is_simd_type<T>::value,int>::type=0>
    inline void subtract(const Matrix<T,Alloc>& a,
                         const Matrix<T,Alloc>& b,
                         Matrix<T,Alloc>& c) {
      assert( (a.rows() == b.rows()) &&
              (a.cols() == b.cols()) );



      if (is_aligned_alloc<Alloc>::value) {        
#ifdef __AVX512F__
        subSIMD<T,Alloc,typename avx512_traits<T>::reg_type>(a,b,c);
#elif  __AVX__
        subSIMD<T,Alloc,typename avx_traits<T>::reg_type>(a,b,c);
#elif  __SSE2__
        subSIMD<T,Alloc,typename sse2_traits<T>::reg_type>(a,b,c);
#else
        ::silent::fallback::subtract(a,b,c);
#endif
      } else { // allocator seems to be unaligned
        ::silent::fallback::subtract(a,b,c);
      }
    }

    // Non-SIMD types such as complex
    template<typename T,
             class Alloc,
             typename std::enable_if<!is_simd_type<T>::value,int>::type = 0>
    inline void subtract(const Matrix<T,Alloc>& a,
                    const Matrix<T,Alloc>& b,
                    Matrix<T,Alloc>& c) {
      
      ::silent::fallback::subtract(a,b,c);
    }

    // In-place implementation a = a-b
    template<typename T,class Alloc>
    inline void subtract(Matrix<T,Alloc>& a,
                         const Matrix<T,Alloc>& b) {

      subtract(a,b,a);
    }


    /*
     * Matrix Multiplication
     */

    // On-copy implementation c=a*b
    template<typename T,class Alloc,typename regType>
    inline void mulSIMD(const Matrix<T,Alloc>& a, 
                        const Matrix<T,Alloc>& b,
                        Matrix<T,Alloc>& c) {

      // This method is instantiated with unaligned allocators.  We
      // allow the instantiation although externally this is never
      // called unaligned
      static_assert(!extract_alignment<Alloc>::aligned ||
        (extract_alignment<Alloc>::value >= sizeof(regType)),
        "Insufficient alignment for the registers used");
      const size_t btentries = b.rows()*b.dcols(); //Takes amount of entries in b
      c.allocate(a.rows(),b.cols()); //Allocates C
      const size_t ctentries = c.rows()*c.dcols();//Takes amoun of entries in c
      regType* here        = reinterpret_cast<regType*>(c.data());//Takes pointer to chunk of C
      regType* cRowStart = here; //Copies the pointer start of C
      const size_t  blocks = ( btentries*sizeof(T) + (sizeof(regType)-1) )/
        sizeof(regType); //Calculates amount of blocks for b
      const size_t  clocks = ( ctentries*sizeof(T) + (sizeof(regType)-1) )/
        sizeof(regType); //Calculates amount of blocks for c
      const regType* bptr  = reinterpret_cast<const regType*>(b.data()); //Takes pointer to chunk of b
      size_t columnBlocks = blocks/b.rows();//Calculates amount of chunks per column for b
      size_t columnClocks = clocks/c.rows();//Calculates amount of chunks per column for c
      const regType* bRowEnd   = bptr + columnBlocks;//Pointer to first chunk of the next row
      const regType* bRowStart   = bptr;

      for(size_t i = 0; i < a.rows(); ++i){//for each row in a
        for(size_t j = 0; j < a.cols(); ++j){ //for each element in that a row
          for (;bptr!=bRowEnd;) {
            T factor = a[i][j];
            if(j == 0){ //Avoid allocation
              *here = scale<T, regType>(factor, *bptr++);//First time just sets the calculation
            }
            else{
              *here = mm_add<T>(*here, scale<T, regType>(factor, *bptr++));//Other times it is added
            }
            here++;//Moves to right in c row
          }
          here = cRowStart;//Restarts position in the c matrix, we need to take precaution on this cause it resets it even when it ends
          bRowEnd += columnBlocks;//moves end to next row of b
        }
        cRowStart += columnClocks; //To calculate next row
        bptr = bRowStart; //reset bptr to start of matrix
        bRowEnd = bRowStart + columnBlocks;//moves next end to the first block of second row
        here = cRowStart; //Sets the pointer for c in the its next row
      }
    }//function mulSIMD

    // In-copy implementation c=a*b
    template<typename T,class Alloc,
       typename std::enable_if<is_simd_type<T>::value,int>::type=0>
    inline void multiply(const Matrix<T,Alloc>& a,
                         const Matrix<T,Alloc>& b,
                         Matrix<T,Alloc>& c) {
      assert(a.cols() == b.rows());
      if (is_aligned_alloc<Alloc>::value && extract_alignment<Alloc>::row_aligned) {      
#ifdef __AVX512F__
        mulSIMD<T,Alloc,typename avx512_traits<T>::reg_type>(a,b,c);
#elif  __AVX__
        mulSIMD<T,Alloc,typename avx_traits<T>::reg_type>(a,b,c);
#elif  __SSE2__
        mulSIMD<T,Alloc,typename sse2_traits<T>::reg_type>(a,b,c);
#else
        ::silent::fallback::multiply(a,b,c);
#endif
      } else { // allocator seems to be unaligned
        ::silent::fallback::multiply(a,b,c);
      }
    }

    // Non-SIMD types such as complex
    template<typename T,
             class Alloc,
             typename std::enable_if<!is_simd_type<T>::value,int>::type = 0>
    inline void multiply(const Matrix<T,Alloc>& a,
                    const Matrix<T,Alloc>& b,
                    Matrix<T,Alloc>& c) {
      ::silent::fallback::multiply(a,b,c);
    }

    // In-place implementation a = a*b
    template<typename T,class Alloc>
    inline void multiply(Matrix<T,Alloc>& a,
                         const Matrix<T,Alloc>& b) {
      Matrix<T,Alloc> c;
      multiply(a,b,c);
      a = c;
    }


    /*
     * Scalar Multiplication
     */

    // On-copy implementation c=a*b
    template<typename T,class Alloc,typename regType>
    inline void mulSIMD(const Matrix<T,Alloc>& a, 
                        const T& b,
                        Matrix<T,Alloc>& c) {

      // This method is instantiated with unaligned allocators.  We
      // allow the instantiation although externally this is never
      // called unaligned
      static_assert(!extract_alignment<Alloc>::aligned ||
        (extract_alignment<Alloc>::value >= sizeof(regType)),
        "Insufficient alignment for the registers used");
      
      const size_t tentries = a.rows()*a.dcols();
      c.allocate(a.rows(),a.cols());
      regType vectorB = mm_set1<T, regType>(b);
      regType* here        = reinterpret_cast<regType*>(c.data());
      const size_t  blocks = ( tentries*sizeof(T) + (sizeof(regType)-1) )/
        sizeof(regType);
      regType *const end   = here + blocks;
      const regType* aptr  = reinterpret_cast<const regType*>(a.data());
      
      for (;here!=end;) {
        *here++ = mm_mul<T, regType>(*aptr++,vectorB);
      }
      
    }//function mulSIMD

    // In-copy implementation c=a*b
    template<typename T,class Alloc,
       typename std::enable_if<is_simd_type<T>::value,int>::type=0>
    inline void multiply(const Matrix<T,Alloc>& a,
                         const T& b,
                         Matrix<T,Alloc>& c) {
      if (is_aligned_alloc<Alloc>::value && extract_alignment<Alloc>::row_aligned) {      
#ifdef __AVX512F__
        mulSIMD<T,Alloc,typename avx512_traits<T>::reg_type>(a,b,c);
#elif  __AVX__
        mulSIMD<T,Alloc,typename avx_traits<T>::reg_type>(a,b,c);
#elif  __SSE2__
        mulSIMD<T,Alloc,typename sse2_traits<T>::reg_type>(a,b,c);
#else
        ::silent::fallback::multiply(a,b,c);
#endif
      } else { // allocator seems to be unaligned
        ::silent::fallback::multiply(a,b,c);
      }
    }

    // Non-SIMD types such as complex
    template<typename T,
             class Alloc,
             typename std::enable_if<!is_simd_type<T>::value,int>::type = 0>
    inline void multiply(const Matrix<T,Alloc>& a,
                    const T& b,
                    Matrix<T,Alloc>& c) {
      ::silent::fallback::multiply(a,b,c);
    }

    // In-place implementation a = a*b
    template<typename T,class Alloc>
    inline void multiply(Matrix<T,Alloc>& a,
                         const T& b) {
      multiply(a,b,a);
    }


  } // namespace simd


  // The arithmetic implementation (aimpl) namespace
  // dispatches to the corresponding methods
#ifdef _ENABLE_SIMD
  namespace aimpl=simd;
#else
  namespace aimpl=fallback;
#endif
  
} // namespace silent

#endif
