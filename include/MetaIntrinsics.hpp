/*
 * Copyright (C) 2017 
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Pablo Alvarado
 * @Date:   28.12.2017
 */

#ifndef SILENT_META_INTRINSICS_HPP
#define SILENT_META_INTRINSICS_HPP

#include "Intrinsics.hpp"

namespace silent {
	namespace simd{

	/*
	 * Vector filling with constant
	*/

    template<typename T,class regType>
    regType mm_set1(T); 
    
#ifdef __AVX512F__
    template<>
    inline __m512d __attribute__((__always_inline__))
    mm_set1<double>(double a) {
      return _mm512_set1_pd(a);
    }
    template<>
    inline __m512 __attribute__((__always_inline__))
    mm_set1<float>(float a) {
      return _mm512_set1_ps(a);
    }
    template<>
    inline __m512i __attribute__((__always_inline__))
    mm_set1<uint64_t>(uint64_t a) {
      return _mm512_set1_epi64(a);
    }
    template<>
    inline __m512i __attribute__((__always_inline__))
    mm_set1<int64_t>(int64_t a) {
      return _mm512_set1_epi64(a);
    }
    template<>
    inline __m512i __attribute__((__always_inline__))
    mm_set1<uint32_t>(uint32_t a) {
      return _mm512_set1_epi32(a);
    }
    template<>
    inline __m512i __attribute__((__always_inline__))
    mm_set1<int32_t>(int32_t a) {
      return _mm512_set1_epi32(a);
    }
    template<>
    inline __m512i __attribute__((__always_inline__))
    mm_set1<uint16_t>(uint16_t a) {
      return _mm512_set1_epi16(a);
    }
    template<>
    inline __m512i __attribute__((__always_inline__))
    mm_set1<int16_t>(int16_t a) {
      return _mm512_set1_epi16(a);
    }
    template<>
    inline __m512i __attribute__((__always_inline__))
    mm_set1<uint8_t>(uint8_t a) {
      return _mm512_set1_epi8(a);
    }
    template<>
    inline __m512i __attribute__((__always_inline__))
    mm_set1<int8_t>(int8_t a) {
      return _mm512_set1_epi8(a);
    }
#elif defined __AVX__
    template<>
    inline __m256d __attribute__((__always_inline__))
    mm_set1<double>(double a) {
      return _mm256_set1_pd(a);
    }
    template<>
    inline __m256 __attribute__((__always_inline__))
    mm_set1<float>(float a) {
      return _mm256_set1_ps(a);
    }
    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_set1<uint64_t>(uint64_t a) {
      return _mm256_set1_epi64x(a);
    }
    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_set1<int64_t>(int64_t a) {
      return _mm256_set1_epi64x(a);
    }
    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_set1<uint32_t>(uint32_t a) {
      return _mm256_set1_epi32(a);
    }
    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_set1<int32_t>(int32_t a) {
      return _mm256_set1_epi32(a);
    }
    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_set1<uint16_t>(uint16_t a) {
      return _mm256_set1_epi16(a);
    }
    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_set1<int16_t>(int16_t a) {
      return _mm256_set1_epi16(a);
    }
    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_set1<uint8_t>(uint8_t a) {
      return _mm256_set1_epi8(a);
    }
    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_set1<int8_t>(int8_t a) {
      return _mm256_set1_epi8(a);
    }
#elif  defined __SSE2__
    template<>
    inline __m128d __attribute__((__always_inline__))
    mm_set1<double>(double a) {
      return _mm128_set1_pd(a);
    }
    template<>
    inline __m128 __attribute__((__always_inline__))
    mm_set1<float>(float a) {
      return _mm128_set1_ps(a);
    }
    template<>
    inline __m128i __attribute__((__always_inline__))
    mm_set1<std::uint64_t>(uint64_t a) {
      return _mm128_set1_epi64(a);
    }
    template<>
    inline __m128i __attribute__((__always_inline__))
    mm_set1<std::int64_t>(int64_t a) {
      return _mm128_set1_epi64(a);
    }
    template<>
    inline __m128i __attribute__((__always_inline__))
    mm_set1<std::uint32_t>(uint32_t a) {
      return _mm128_set1_epi32(a);
    }
    template<>
    inline __m128i __attribute__((__always_inline__))
    mm_set1<std::int32_t>(int32_t a) {
      return _mm128_set1_epi32(a);
    }
    template<>
    inline __m128i __attribute__((__always_inline__))
    mm_set1<std::uint16_t>(uint16_t a) {
      return _mm128_set1_epi16(a);
    }
    template<>
    inline __m128i __attribute__((__always_inline__))
    mm_set1<std::int16_t>(int16_t a) {
      return _mm128_set1_epi16(a);
    }
    template<>
    inline __m128i __attribute__((__always_inline__))
    mm_set1<std::uint8_t>(uint8_t a) {
      return _mm128_set1_epi8(a);
    }
    template<>
    inline __m128i __attribute__((__always_inline__))
    mm_set1<std::int8_t>(int8_t a) {
      return _mm128_set1_epi8(a);
    }
#endif

/*
* Multiplication for each element
*/

    template<typename T,class regType>
    regType mm_mul(regType, regType); 

#ifdef __AVX512F__
    template<>
    inline __m512d __attribute__((__always_inline__))
    mm_mul<double>(__m512d a,__m512d b) {
      return _mm512_mul_pd(a,b);
    }
    template<>
    inline __m512 __attribute__((__always_inline__))
    mm_mul<float>(__m512 a,__m512 b) {
      return _mm512_mul_ps(a,b);
    }
    template<>
    inline __m512i __attribute__((__always_inline__))
    mm_mul<uint64_t>(__m512i a,__m512i b) {
      return _mm512_mul_epi64(a,b);
    }
    template<>
    inline __m512i __attribute__((__always_inline__))
    mm_mul<int64_t>(__m512i a,__m512i b) {
      return _mm512_mul_epi64(a,b);
    }
    template<>
    inline __m512i __attribute__((__always_inline__))
    mm_mul<uint32_t>(__m512i a,__m512i b) {
      return _mm512_mul_epi32(a,b);
    }
    template<>
    inline __m512i __attribute__((__always_inline__))
    mm_mul<int32_t>(__m512i a,__m512i b) {
      return _mm512_mul_epi32(a,b);
    }
    template<>
    inline __m512i __attribute__((__always_inline__))
    mm_mul<uint16_t>(__m512i a,__m512i b) {
      return _mm512_mul_epi16(a,b);
    }
    template<>
    inline __m512i __attribute__((__always_inline__))
    mm_mul<int16_t>(__m512i a,__m512i b) {
      return _mm512_mul_epi16(a,b);
    }
    template<>
    inline __m512i __attribute__((__always_inline__))
    mm_mul<uint8_t>(__m512i a,__m512i b) {
      return _mm512_mul_epi8(a,b);
    }
    template<>
    inline __m512i __attribute__((__always_inline__))
    mm_mul<int8_t>(__m512i a,__m512i b) {
      return _mm512_mul_epi8(a,b);
    }
#elif defined __AVX__
    template<>
    inline __m256d __attribute__((__always_inline__))
    mm_mul<double>(__m256d a,__m256d b) {
      return _mm256_mul_pd(a,b);
    }
    template<>
    inline __m256 __attribute__((__always_inline__))
    mm_mul<float>(__m256 a,__m256 b) {
      return _mm256_mul_ps(a,b);
    }
    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_mul<uint64_t>(__m256i a,__m256i b) {
      return _mm256_mullo_epi64(a,b);
    }
    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_mul<int64_t>(__m256i a,__m256i b) {
      return _mm256_mullo_epi64(a,b);
    }
    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_mul<uint32_t>(__m256i a,__m256i b) {
      return _mm256_mullo_epi32(a,b);
    }
    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_mul<int32_t>(__m256i a,__m256i b) {
      return _mm256_mullo_epi32(a,b);
    }
    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_mul<uint16_t>(__m256i a,__m256i b) {
      return _mm256_mullo_epi16(a,b);
    }
    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_mul<int16_t>(__m256i a,__m256i b) {
      return _mm256_mullo_epi16(a,b);
    }
    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_mul<uint8_t>(__m256i a,__m256i b) {
      return _mm256_mullo_epi16(a,b); //High risk of nor working
    }
    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_mul<int8_t>(__m256i a,__m256i b) {
      return _mm256_mullo_epi16(a,b); //High risk of nor working
    }
#elif  defined __SSE2__
    template<>
    inline __m128d __attribute__((__always_inline__))
    mm_mul<double>(__m128d a,__m128d b) {
      return _mm_mul_pd(a,b);
    }
    template<>
    inline __m128 __attribute__((__always_inline__))
    mm_mul<float>(__m128 a,__m128 b) {
      return _mm_mul_ps(a,b);
    }
    template<>
    inline __m128i __attribute__((__always_inline__))
    mm_mul<std::uint64_t>(__m128i a,__m128i b) {
      return _mm_mul_epi64(a,b);
    }
    template<>
    inline __m128i __attribute__((__always_inline__))
    mm_mul<std::int64_t>(__m128i a,__m128i b) {
      return _mm_mul_epi64(a,b);
    }
    template<>
    inline __m128i __attribute__((__always_inline__))
    mm_mul<std::uint32_t>(__m128i a,__m128i b) {
      return _mm_mul_epi32(a,b);
    }
    template<>
    inline __m128i __attribute__((__always_inline__))
    mm_mul<std::int32_t>(__m128i a,__m128i b) {
      return _mm_mul_epi32(a,b);
    }
    template<>
    inline __m128i __attribute__((__always_inline__))
    mm_mul<std::uint16_t>(__m128i a,__m128i b) {
      return _mm_mul_epi16(a,b);
    }
    template<>
    inline __m128i __attribute__((__always_inline__))
    mm_mul<std::int16_t>(__m128i a,__m128i b) {
      return _mm_mul_epi16(a,b);
    }
    template<>
    inline __m128i __attribute__((__always_inline__))
    mm_mul<std::uint8_t>(__m128i a,__m128i b) {
      return _mm_mul_epi8(a,b);
    }
    template<>
    inline __m128i __attribute__((__always_inline__))
    mm_mul<std::int8_t>(__m128i a,__m128i b) {
      return _mm_mul_epi8(a,b);
    }
#endif

    /*
     * Subtraction
     */


    /// We wrap the intrinsics methods to be polymorphic versions
    template<typename T,class regType>
    regType mm_sub(regType,regType); // We don't implement this to cause, at
                                     // least, a linker error if this version is
                                     // used.
    //{
    // Generic function should never be called.
    // If it is called, then some SIMD chaos is going on...
    
    // A way to cause a compile time error would be better
    // throw std::bad_function_call();
    // return regType();
    //}

#ifdef __AVX512F__
    template<>
    inline __m512d __attribute__((__always_inline__))
    mm_sub<double>(__m512d a,__m512d b) {
      return _mm512_sub_pd(a,b);
    }
    template<>
    inline __m512 __attribute__((__always_inline__))
    mm_sub<float>(__m512 a,__m512 b) {
      return _mm512_sub_ps(a,b);
    }
    template<>
    inline __m512i __attribute__((__always_inline__))
    mm_sub<uint64_t>(__m512i a,__m512i b) {
      return _mm512_sub_epi64(a,b);
    }
    template<>
    inline __m512i __attribute__((__always_inline__))
    mm_sub<int64_t>(__m512i a,__m512i b) {
      return _mm512_sub_epi64(a,b);
    }
    template<>
    inline __m512i __attribute__((__always_inline__))
    mm_sub<uint32_t>(__m512i a,__m512i b) {
      return _mm512_sub_epi32(a,b);
    }
    template<>
    inline __m512i __attribute__((__always_inline__))
    mm_sub<int32_t>(__m512i a,__m512i b) {
      return _mm512_sub_epi32(a,b);
    }
    template<>
    inline __m512i __attribute__((__always_inline__))
    mm_sub<uint16_t>(__m512i a,__m512i b) {
      return _mm512_sub_epi16(a,b);
    }
    template<>
    inline __m512i __attribute__((__always_inline__))
    mm_sub<int16_t>(__m512i a,__m512i b) {
      return _mm512_sub_epi16(a,b);
    }
    template<>
    inline __m512i __attribute__((__always_inline__))
    mm_sub<uint8_t>(__m512i a,__m512i b) {
      return _mm512_sub_epi8(a,b);
    }
    template<>
    inline __m512i __attribute__((__always_inline__))
    mm_sub<int8_t>(__m512i a,__m512i b) {
      return _mm512_sub_epi8(a,b);
    }
#elif defined __AVX__
    template<>
    inline __m256d __attribute__((__always_inline__))
    mm_sub<double>(__m256d a,__m256d b) {
      return _mm256_sub_pd(a,b);
    }
    template<>
    inline __m256 __attribute__((__always_inline__))
    mm_sub<float>(__m256 a,__m256 b) {
      return _mm256_sub_ps(a,b);
    }
    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_sub<uint64_t>(__m256i a,__m256i b) {
      return _mm256_sub_epi64(a,b);
    }
    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_sub<int64_t>(__m256i a,__m256i b) {
      return _mm256_sub_epi64(a,b);
    }
    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_sub<uint32_t>(__m256i a,__m256i b) {
      return _mm256_sub_epi32(a,b);
    }
    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_sub<int32_t>(__m256i a,__m256i b) {
      return _mm256_sub_epi32(a,b);
    }
    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_sub<uint16_t>(__m256i a,__m256i b) {
      return _mm256_sub_epi16(a,b);
    }
    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_sub<int16_t>(__m256i a,__m256i b) {
      return _mm256_sub_epi16(a,b);
    }
    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_sub<uint8_t>(__m256i a,__m256i b) {
      return _mm256_sub_epi8(a,b);
    }
    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_sub<int8_t>(__m256i a,__m256i b) {
      return _mm256_sub_epi8(a,b);
    }
#elif defined __SSE2__
    template<>
    inline __m128d __attribute__((__always_inline__))
    mm_sub<double>(__m128d a,__m128d b) {
      return _mm_sub_pd(a,b);
    }
    template<>
    inline __m128 __attribute__((__always_inline__))
    mm_sub<float>(__m128 a,__m128 b) {
      return _mm_sub_ps(a,b);
    }
    template<>
    inline __m128i __attribute__((__always_inline__))
    mm_sub<std::uint64_t>(__m128i a,__m128i b) {
      return _mm_sub_epi64(a,b);
    }
    template<>
    inline __m128i __attribute__((__always_inline__))
    mm_sub<std::int64_t>(__m128i a,__m128i b) {
      return _mm_sub_epi64(a,b);
    }
    template<>
    inline __m128i __attribute__((__always_inline__))
    mm_sub<std::uint32_t>(__m128i a,__m128i b) {
      return _mm_sub_epi32(a,b);
    }
    template<>
    inline __m128i __attribute__((__always_inline__))
    mm_sub<std::int32_t>(__m128i a,__m128i b) {

      return _mm_sub_epi32(a,b);
    }
    template<>
    inline __m128i __attribute__((__always_inline__))
    mm_sub<std::uint16_t>(__m128i a,__m128i b) {
      return _mm_sub_epi16(a,b);
    }
    template<>
    inline __m128i __attribute__((__always_inline__))
    mm_sub<std::int16_t>(__m128i a,__m128i b) {
      return _mm_sub_epi16(a,b);
    }
    template<>
    inline __m128i __attribute__((__always_inline__))
    mm_sub<std::uint8_t>(__m128i a,__m128i b) {
      return _mm_sub_epi8(a,b);
    }
    template<>
    inline __m128i __attribute__((__always_inline__))
    mm_sub<std::int8_t>(__m128i a,__m128i b) {
      return _mm_sub_epi8(a,b);
    }
#endif

/*
     * Sum
     */

    /*
     * The following code exemplifies how to manually accelerate code using
     * SIMD instructions.  However, for the simple element-wise algorithms
     * like sum or subtraction, modern compilers can automatically vectorize
     * the code, as the benchmarks show.
     */


    /// We wrap the intrinsics methods to be polymorphic versions
    template<typename T,class regType>
    regType mm_add(regType,regType); // We don't implement this to cause, at
                                     // least, a linker error if this version is
                                     // used.
    //{
    // Generic function should never be called.
    // If it is called, then some SIMD chaos is going on...
    
    // A way to cause a compile time error would be better
    // throw std::bad_function_call();
    // return regType();
    //}

    
#ifdef __AVX512F__
    template<>
    inline __m512d __attribute__((__always_inline__))
    mm_add<double>(__m512d a,__m512d b) {
      return _mm512_add_pd(a,b);
    }
    template<>
    inline __m512 __attribute__((__always_inline__))
    mm_add<float>(__m512 a,__m512 b) {
      return _mm512_add_ps(a,b);
    }
    template<>
    inline __m512i __attribute__((__always_inline__))
    mm_add<uint64_t>(__m512i a,__m512i b) {
      return _mm512_add_epi64(a,b);
    }
    template<>
    inline __m512i __attribute__((__always_inline__))
    mm_add<int64_t>(__m512i a,__m512i b) {
      return _mm512_add_epi64(a,b);
    }
    template<>
    inline __m512i __attribute__((__always_inline__))
    mm_add<uint32_t>(__m512i a,__m512i b) {
      return _mm512_add_epi32(a,b);
    }
    template<>
    inline __m512i __attribute__((__always_inline__))
    mm_add<int32_t>(__m512i a,__m512i b) {
      return _mm512_add_epi32(a,b);
    }
    template<>
    inline __m512i __attribute__((__always_inline__))
    mm_add<uint16_t>(__m512i a,__m512i b) {
      return _mm512_add_epi16(a,b);
    }
    template<>
    inline __m512i __attribute__((__always_inline__))
    mm_add<int16_t>(__m512i a,__m512i b) {
      return _mm512_add_epi16(a,b);
    }
    template<>
    inline __m512i __attribute__((__always_inline__))
    mm_add<uint8_t>(__m512i a,__m512i b) {
      return _mm512_add_epi8(a,b);
    }
    template<>
    inline __m512i __attribute__((__always_inline__))
    mm_add<int8_t>(__m512i a,__m512i b) {
      return _mm512_add_epi8(a,b);
    }
#elif defined __AVX__
    template<>
    inline __m256d __attribute__((__always_inline__))
    mm_add<double>(__m256d a,__m256d b) {
      return _mm256_add_pd(a,b);
    }
    template<>
    inline __m256 __attribute__((__always_inline__))
    mm_add<float>(__m256 a,__m256 b) {
      return _mm256_add_ps(a,b);
    }
    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_add<uint64_t>(__m256i a,__m256i b) {
      return _mm256_add_epi64(a,b);
    }
    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_add<int64_t>(__m256i a,__m256i b) {
      return _mm256_add_epi64(a,b);
    }
    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_add<uint32_t>(__m256i a,__m256i b) {
      return _mm256_add_epi32(a,b);
    }
    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_add<int32_t>(__m256i a,__m256i b) {
      return _mm256_add_epi32(a,b);
    }
    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_add<uint16_t>(__m256i a,__m256i b) {
      return _mm256_add_epi16(a,b);
    }
    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_add<int16_t>(__m256i a,__m256i b) {
      return _mm256_add_epi16(a,b);
    }
    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_add<uint8_t>(__m256i a,__m256i b) {
      return _mm256_add_epi8(a,b);
    }
    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_add<int8_t>(__m256i a,__m256i b) {
      return _mm256_add_epi8(a,b);
    }
#elif  defined __SSE2__
    template<>
    inline __m128d __attribute__((__always_inline__))
    mm_add<double>(__m128d a,__m128d b) {
      return _mm_add_pd(a,b);
    }
    template<>
    inline __m128 __attribute__((__always_inline__))
    mm_add<float>(__m128 a,__m128 b) {
      return _mm_add_ps(a,b);
    }
    template<>
    inline __m128i __attribute__((__always_inline__))
    mm_add<std::uint64_t>(__m128i a,__m128i b) {
      return _mm_add_epi64(a,b);
    }
    template<>
    inline __m128i __attribute__((__always_inline__))
    mm_add<std::int64_t>(__m128i a,__m128i b) {
      return _mm_add_epi64(a,b);
    }
    template<>
    inline __m128i __attribute__((__always_inline__))
    mm_add<std::uint32_t>(__m128i a,__m128i b) {
      return _mm_add_epi32(a,b);
    }
    template<>
    inline __m128i __attribute__((__always_inline__))
    mm_add<std::int32_t>(__m128i a,__m128i b) {
      return _mm_add_epi32(a,b);
    }
    template<>
    inline __m128i __attribute__((__always_inline__))
    mm_add<std::uint16_t>(__m128i a,__m128i b) {
      return _mm_add_epi16(a,b);
    }
    template<>
    inline __m128i __attribute__((__always_inline__))
    mm_add<std::int16_t>(__m128i a,__m128i b) {
      return _mm_add_epi16(a,b);
    }
    template<>
    inline __m128i __attribute__((__always_inline__))
    mm_add<std::uint8_t>(__m128i a,__m128i b) {
      return _mm_add_epi8(a,b);
    }
    template<>
    inline __m128i __attribute__((__always_inline__))
    mm_add<std::int8_t>(__m128i a,__m128i b) {
      return _mm_add_epi8(a,b);
    }
#endif
	}//namespace simd
}//namespace silent
#endif