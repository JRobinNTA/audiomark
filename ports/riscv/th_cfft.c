#include "th_cfft.h"

/* The types here are kept explicit for readability */
#if defined(__riscv_p) && !defined(__OPTIMIZE__)

static void riscv_cfft_radix4by2_inverse(riscv_cfft_instance *p_instance, q31_t *p_Src, uint32_t fftLen){

    uint32_t     n2;
    q31_t       *pIn0;
    q31_t       *pIn1;
    const q31_t *pCoef = p_instance->pTwiddle;
    uint32_t     blkCnt;
    q31x2_t    vecIn0, vecIn1, vecSum, vecDiff, vecCmplxTmp;
    q31x2_t    vecTmp, vecTmpx, vecTw, vecSwp, vecSwpx;

    n2 = fftLen >> 1;
    pIn0 = pSrc;
    pIn1 = pSrc + fftLen;

    blkCnt = n2;

    while (blkCnt > 0U)
    {
        vecIn0 = __riscv_pload_i32x2(pIn0);
        vecIn1 = __riscv_pload_i32x2(pIn1);

        vecIn0 = __riscv_psra_s_i32x2(vecIn0, 1);
        vecIn1 = __riscv_psra_s_i32x2(vecIn1, 1);
        vecSum = __riscv_paadd_i32x2(vecIn0, vecIn1);
        __riscv_pstore_i32x2(pIn0, vecSum);
        pIn0 += 2;

        vecTw = __riscv_pload_i32x2(pCoef);
        pCoef += 2;
        vecDiff = __riscv_pasub_i32x2(vecIn0, vecIn1);

        /* { vecDiff[0]*vecTw[0], vecDiff[1]*vecTw[1] } */
        vecTmp = __riscv_pmulqr_i32x2(vecDiff,vecTw);
        /* { vecDiff[0]*vecTw[1], vecDiff[1]*vecTw[0] } */
        vecTmpx = __riscv_pmulqr_i32x2(vecDiff, __riscv_ppairoe_i32x2(vecTw, vecTw));
        /* vecSwp  = { vecDiff[0]*vecTw[0], vecDiff[1]*vecTw[0] } */
        vecSwp = __riscv_ppaireo_i32x2(vecTmp, vecTmpx);
        /* vecSwpx = { vecDiff[1]*vecTw[1], vecDiff[0]*vecTw[1] } */
        vecSwpx = __riscv_ppairoe_i32x2(vecTmp, vecTmpx);
        /* vecCmplxTmp = { vecDiff[0]*vecTw[0] - vecDiff[1]*vecTw[1],
           vecDiff[1]*vecTw[0] + vecDiff[0]*vecTw[1]} */
        vecCmplxTmp = __riscv_psa_x_i32x2(vecSwp, vecSwpx);
        __riscv_pstore_i32x2(pIn1, vecCmplxTmp);
        pIn1 += 2;

        blkCnt--;
    }

    riscv_radix4_butterfly_q31(p_instance, pSrc, n2);

    riscv_radix4_butterfly_q31(p_instance, pSrc + fftLen, n2);

    pIn0 = pSrc;
    blkCnt = (fftLen << 1) >> 2;
    while (blkCnt > 0U)
    {
        vecIn0 = __riscv_pload_i32x2(pIn0);
        vecIn0 = __riscv_psll_s_i32x2(vecIn0, 1);
        __riscv_pstore_i32x2(pIn0, vecIn0);
        pIn0 += 2;
        blkCnt--;
    }
    /* No tail handling required since fftLen is always a multiple of 2 */
}

static void riscv_cfft_radix4by2(riscv_cfft_instance *p_instance, q31_t *p_Src, uint32_t fftLen){

}

static void riscv_bitreversal_32_inpl(uint32_t *p_Src, uint16_t bitRevLength, const uint16_t *pBitRevTable){

}

#else

static void riscv_cfft_radix4by2_inverse(riscv_cfft_instance *p_instance, float32_t *p_Src, uint32_t fftLen){

}

static void riscv_cfft_radix4by2(riscv_cfft_instance *p_instance, float32_t *p_Src, uint32_t fftLen){

}

static void riscv_bitreversal_32_inpl(uint32_t *p_Src, uint16_t bitRevLength, const uint16_t *pBitRevTable){

}
#endif
