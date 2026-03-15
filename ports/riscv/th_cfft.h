#ifndef __TH_CFFT_H
#define __TH_CFFT_H

#include "th_types.h"


#if defined(__riscv_p)

#define RISCVBITREVINDEXTABLE_128_TABLE_LENGTH ((uint16_t)112)

/* vecC = { vecA[low]*vecB[low] - vecA[high]*vecB[high],
           vecA[high]*vecB[low] + vecA[low]*vecB[high]} */
#define RISCV_PEXT_CMPLX_MULT_FX_AxB(_vecA, _vecB, _vecC) \
    do { \
        q31x2_t _prod1 = __riscv_pmulqr_i32x2(_vecA, _vecB); \
        q31x2_t _swapB = __riscv_ppairoe_i32x2(_vecB, _vecB); \
        q31x2_t _prod2 = __riscv_pmulqr_i32x2(_vecA, _swapB); \
        _vecC = __riscv_psa_x_i32x2( \
            __riscv_ppaireo_i32x2(_prod1, _prod2), \
            __riscv_ppairoe_i32x2(_prod1, _prod2)  \
        ); \
    } while (0)

/* vecC = { vecA[low]*vecB[low] + vecA[high]*vecB[high],
           vecA[high]*vecB[low] - vecA[low]*vecB[high] } */
#define RISCV_PEXT_CMPLX_MULT_FX_AxConjB(_vecA, _vecB, _vecC) \
    do { \
        q31x2_t _prod1 = __riscv_pmulqr_i32x2(_vecA, _vecB); \
        q31x2_t _swapB = __riscv_ppairoe_i32x2(_vecB, _vecB); \
        q31x2_t _prod2 = __riscv_pmulqr_i32x2(_vecA, _swapB); \
        _vecC = __riscv_pas_x_i32x2( \
            __riscv_ppaireo_i32x2(_prod1, _prod2), \
            __riscv_ppairoe_i32x2(_prod1, _prod2)  \
        ); \
    } while (0)

/* vecC = { vecA[low] - vecB[high],
           vecA[high] + vecB[low] } */
#define RISCV_PEXT_CMPLX_ADD_FX_A_ixB(_vecA, _vecB, _vecC) \
    do { \
        _vecC = __riscv_paas_x_i32x2(_vecA, _vecB); \
    } while (0)

/* vecC = { vecA[low] + vecB[high],
           vecA[high] - vecB[low] } */
#define RISCV_PEXT_CMPLX_SUB_FX_A_ixB(_vecA, _vecB, _vecC) \
    do { \
        _vecC = __riscv_pasa_x_i32x2(_vecA, _vecB); \
    } while (0)

extern const dsp_scalar_t twiddleCoef_128[192];

extern const uint32_t rearranged_twiddle_tab_stride1_arr_64_q31[3];

extern const uint32_t rearranged_twiddle_tab_stride2_arr_64_q31[3];

extern const uint32_t rearranged_twiddle_tab_stride3_arr_64_q31[3];

extern const q31_t rearranged_twiddle_stride1_64_q31[40];

extern const q31_t rearranged_twiddle_stride2_64_q31[40];

extern const q31_t rearranged_twiddle_stride3_64_q31[40];

void riscv_cfft_radix4by2_inverse(riscv_cfft_instance *p_instance, q31_t *pSrc, uint32_t fftLen);

void riscv_cfft_radix4by2(riscv_cfft_instance *p_instance, q31_t *pSrc, uint32_t fftLen);

void riscv_radix4_butterfly_inverse_q31(const riscv_cfft_instance *S, q31_t   *pSrc, uint32_t fftLen);

void riscv_radix4_butterfly_q31(const riscv_cfft_instance *S, q31_t   *pSrc, uint32_t fftLen);

void riscv_bitreversal_32_inpl(uint32_t* p_Src, uint16_t bitRevLength, const uint16_t *pBitRevTable);

void riscv_cfft_init(riscv_cfft_instance *p_instance);
#else

#define RISCVBITREVINDEXTABLE_128_TABLE_LENGTH ((uint16_t)208)
extern const dsp_scalar_t twiddleCoef_128[256];
static void riscv_cfft_radix4by2_inverse(riscv_cfft_instance *p_instance, dsp_scalar_t *pSrc, uint32_t fftLen);

static void riscv_cfft_radix4by2(riscv_cfft_instance *p_instance, dsp_scalar_t *pSrc, uint32_t fftLen);

static void riscv_bitreversal_32_inpl(uint32_t* pSrc, uint16_t bitRevLength, const uint16_t *pBitRevTable);

#endif

extern const uint16_t riscvBitRevIndexTable128[RISCVBITREVINDEXTABLE_128_TABLE_LENGTH];

#endif
