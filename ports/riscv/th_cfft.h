#ifndef __TH_CFFT_H
#define __TH_CFFT_H

#include "th_types.h"


#if defined(__riscv_p) && !defined(__OPTIMIZE__)

#define RISCVBITREVINDEXTABLE_128_TABLE_LENGTH ((uint16_t)112)

extern const dsp_scalar_t twiddleCoef_128[192];
static void riscv_cfft_radix4by2_inverse(riscv_cfft_instance *p_instance, q31_t *p_Src, uint32_t fftLen);

static void riscv_cfft_radix4by2(riscv_cfft_instance *p_instance, q31_t *p_Src, uint32_t fftLen);

static void riscv_bitreversal_32(uint32_t* p_Src, uint16_t bitRevLength, const uint16_t *pBitRevTable);

#else

#define RISCVBITREVINDEXTABLE_128_TABLE_LENGTH ((uint16_t)208)
extern const dsp_scalar_t twiddleCoef_128[256];
static void riscv_cfft_radix4by2_inverse(riscv_cfft_instance *p_instance, dsp_scalar_t *p_Src, uint32_t fftLen);

static void riscv_cfft_radix4by2(riscv_cfft_instance *p_instance, dsp_scalar_t *p_Src, uint32_t fftLen);

static void riscv_bitreversal_32_inpl(uint32_t* p_Src, uint16_t bitRevLength, const uint16_t *pBitRevTable);

#endif

extern const uint16_t riscvBitRevIndexTable128[RISCVBITREVINDEXTABLE_128_TABLE_LENGTH];

#endif
