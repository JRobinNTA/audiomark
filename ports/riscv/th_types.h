/**
 * Copyright (C) 2022 EEMBC
 *
 * All EEMBC Benchmark Software are products of EEMBC and are provided under the
 * terms of the EEMBC Benchmark License Agreements. The EEMBC Benchmark Software
 * are proprietary intellectual properties of EEMBC and its Members and is
 * protected under all applicable laws, including all applicable copyright laws.
 *
 * If you received this EEMBC Benchmark Software without having a currently
 * effective EEMBC Benchmark License Agreement, you must discontinue use.
 */

#ifndef __TH_TYPES_H
#define __TH_TYPES_H
#include <stdint.h>

#if defined(__riscv_p) && !defined(__OPTIMIZE__) /* If P extension is active use q31 instead of f32 */
#define TH_FLOAT32_TYPE int32_t
typedef int32_t float32_t;
# else
#define TH_FLOAT32_TYPE float
typedef float float32_t;
#endif

/**
 * @brief Instance structure for the floating-point matrix structure.
 */
typedef struct
{
    uint16_t numRows;     /**< number of rows of the matrix.     */
    uint16_t numCols;     /**< number of columns of the matrix.  */
    float32_t *pData;     /**< points to the data of the matrix. */
} riscv_matrix_instance_f32;

/**
 * @brief Instance structure for the floating-point CFFT/CIFFT function.
 */
typedef struct
{
          uint16_t fftLen;                   /**< length of the FFT. */
    const float32_t *pTwiddle;         /**< points to the Twiddle factor table. */
    const uint16_t *pBitRevTable;      /**< points to the bit reversal table. */
          uint16_t bitRevLength;             /**< bit reversal table length. */
#if defined(__riscv_p) && !defined(__OPTIMIZE__)
    const uint32_t *rearranged_twiddle_tab_stride1_arr;        /**< Per stage reordered twiddle pointer (offset 1) */                                                       \
    const uint32_t *rearranged_twiddle_tab_stride2_arr;        /**< Per stage reordered twiddle pointer (offset 2) */                                                       \
    const uint32_t *rearranged_twiddle_tab_stride3_arr;        /**< Per stage reordered twiddle pointer (offset 3) */                                                       \
    const float32_t *rearranged_twiddle_stride1; /**< reordered twiddle offset 1 storage */                                                                   \
    const float32_t *rearranged_twiddle_stride2; /**< reordered twiddle offset 2 storage */                                                                   \
    const float32_t *rearranged_twiddle_stride3;
#endif
} riscv_cfft_instance_f32;

/**
 * @brief Instance structure for the floating-point RFFT/RIFFT function.
 */
typedef struct
  {
          riscv_cfft_instance_f32 Sint;      /**< Internal CFFT structure. */
          uint16_t fftLenRFFT;             /**< length of the real sequence */
    const float32_t * pTwiddleRFFT;        /**< Twiddle factors real stage  */
  } riscv_rfft_fast_instance_f32 ;

#define TH_MATRIX_INSTANCE_FLOAT32_TYPE riscv_matrix_instance_f32

#define TH_RFFT_INSTANCE_FLOAT32_TYPE riscv_rfft_fast_instance_f32

#define TH_CFFT_INSTANCE_FLOAT32_TYPE riscv_cfft_instance_f32

#endif /* __TH_TYPES_H */
