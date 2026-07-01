/**
 * \file csc.c
 * Analysis of clock skew compensation (CSC) algorithms.
 *
 * \author Kyeong Soo (Joseph) Kim <kyeongsoo.kim@gmail.com>
 *
 * \note The following options are controlled by macro definitions:
 * - CSC_INT_SIZE: The number of bytes for 'i', 'D', and 'A' (4 or 8).
 * - CSC_NO_ITER_COUNT: Turn off iteration couting in iterative algorithms.
 * - CSC_NO_DIV_CHECK: Turn off checking the value of 'A' in division algorithms.
 *
 * \remarks The results are published in the following paper:
 * - Kyeong Soo Kim, "Space-time trade-off in integer linear scaling
 *   rounded to the nearest integer through multiplicative and additive
 *   decomposition," arXiv e-prints arXiv:2605.21400v [cs.DS], May 2026.
 *   [Online]. Available: https://arxiv.org/abs/2605.21400
 */

#include <assert.h>
#include <inttypes.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include "csc.h"

#if defined(__x86_64__) || defined(_M_X64) || defined(__amd64__)
#include <gmp.h>
#include <mpfr.h>
#define IEEE754_B128_PREC 113 // for IEEE 754-2008 binary128 format
#define IEEE754_B256_PREC 237 // for IEEE 754-2008 binary256 format
#define IEEE754_B512_PREC 489 // for IEEE 754-2008 binary512 format

/**
 * \brief CSC based on mpfr FP division.
 */
csc_int_t csc_mpfr(const csc_int_t i, const csc_int_t D, const csc_int_t A, int precision)
{
#ifndef CSC_NO_DIV_CHECK
    if (A == 0) {
        return 0;
    }
#endif
    mpfr_t mpfr_D, mpfr_A;

    mpfr_set_default_prec(precision);
    mpfr_init(mpfr_D);
    mpfr_init(mpfr_A);

     // integers are first converted to double to deal with very large values (e.g., > INT64_MAX)
    mpfr_set_d(mpfr_D, (double)D, MPFR_RNDN);
    mpfr_set_d(mpfr_A, (double)A, MPFR_RNDN);
    mpfr_mul_d(mpfr_D, mpfr_D, (double)i, MPFR_RNDN);
    mpfr_div(mpfr_D, mpfr_D, mpfr_A, MPFR_RNDN);
    mpfr_add_d(mpfr_D, mpfr_D, 0.5, MPFR_RNDN);
    mpfr_floor(mpfr_D, mpfr_D);
    #if CSC_INT_SIZE == 16 && defined(__SIZEOF_INT128__)
    mpz_t mpz_result;
    mpz_init(mpz_result);
    mpfr_get_z(mpz_result, mpfr_D, MPFR_RNDN);
    csc_int_t csc_result;
    mpz_export(&csc_result, NULL, -1, sizeof(csc_result), 0, 0, mpz_result);
    mpz_clear(mpz_result);
    return csc_result;
    #else
    return (csc_int_t) mpfr_get_sj(mpfr_D, MPFR_RNDN);
    #endif
}

/**
 * \brief A wrapper function for CSC based on binary128 FP division.
 */
csc_int_t csc_b128(const csc_int_t i, const csc_int_t D, const csc_int_t A)
{
    return csc_mpfr(i, D, A, IEEE754_B128_PREC);
}

/**
 * \brief A wrapper function for CSC based on binary256 FP division.
 */
csc_int_t csc_b256(const csc_int_t i, const csc_int_t D, const csc_int_t A)
{
    return csc_mpfr(i, D, A, IEEE754_B256_PREC);
}

/**
 * \brief A wrapper function for CSC based on binary512 FP division.
 */
csc_int_t csc_b512(const csc_int_t i, const csc_int_t D, const csc_int_t A)
{
    return csc_mpfr(i, D, A, IEEE754_B512_PREC);
}
#endif // for "#if defined(__x86_64__) || defined(_M_X64) || defined(__amd64__)"

/**
 * \brief CSC based on double-precision FP division.
 */
csc_int_t csc_dp(const csc_int_t i, const csc_int_t D, const csc_int_t A)
{
#ifndef CSC_NO_DIV_CHECK
    if (A == 0) {
        return 0;
    }
#endif
    // return (csc_int_t) floor((i*(double)D/(double)A) + 0.5); // floor() not working for uint64_t on TelosB platform
    return (csc_int_t)((i * (double)D / (double)A) + 0.5);
}

/**
 * \brief CSC based on single-precision FP division.
 */
csc_int_t csc_sp(const csc_int_t i, const csc_int_t D, const csc_int_t A)
{
#ifndef CSC_NO_DIV_CHECK
    if (A == 0) {
        return 0;
    }
#endif
    // return (csc_int_t) floor((i*(float)D/(float)A) + 0.5); // floor() not working for uint64_t on TelosB platform
    return (csc_int_t)((i * (float)D / (float)A) + 0.5);
}

/**
 * \brief CSC based on the extended Bresenham's algorithm with approximate bounds.
 *
 * \remarks For details, refer to the following paper:
 * - K. S. Kim and S. Kang, "Clock skew compensation algorithm immune to
 *   floating-point precision loss," IEEE Commun. Lett., vol. 26, no. 4,
 *   pp. 902-–906, Apr. 2022.
 */
csc_int_t csc_eb1(const csc_int_t i, const csc_int_t D, const csc_int_t A, float epsilon, uint32_t *p_num_iter)
{
#ifndef CSC_NO_DIV_CHECK
    if (A == 0) {
        return 0;
    }
#endif
    csc_int_t da = A;
    csc_int_t db = (D > A) ? (D - A) : D; // Delta a and Delta b
    csc_int_t lb = (csc_int_t)ceilf(i * ((float)db / (float)da - epsilon) - 1);
    csc_int_t ub = (csc_int_t)floorf(i * ((float)db / (float)da + epsilon) + 1);

    // starting point for the iteration
    csc_int_t sp = i - (ub - lb); // x coordinate
    csc_int_t y = lb; // y coordinate
    for (csc_int_t x = sp; x < i; x++) {
        // process overlined triangle down in Eq. (9) of the reference
        y += (2 * (x * db - y * da) >= 0);
    }
    *p_num_iter = ub - lb;
    return (D > A) ? (y + i) : y;
}

/**
 * \brief CSC based on the extended Bresenham's algorithm with theoretical bounds.
 *
 * \remarks For details, refer to the following paper:
 * - S. Kang and K. S. Kim, "Theoretical and practical bounds on the initial value
 *   of clock skew compensation algorithm immune to floating-point precision loss
 *   for resource-constrained wireless sensor nodes," Fiber and Integrated Optics,
 *   vol. 43, no. 3, pp. 111–121, Jun. 24, 2024.
 */
csc_int_t csc_eb2(const csc_int_t i, const csc_int_t D, const csc_int_t A, float epsilon, uint32_t *p_num_iter)
{
#ifndef CSC_NO_DIV_CHECK
    if (A == 0) {
        return 0;
    }
#endif
    csc_int_t da = A;
    csc_int_t db = (D > A) ? (D - A) : D; // Delta a and Delta b

    float u = powf(2, -24);
    float comp1 = 1 + u;
    float comp2 = 1 + 2 * u;
    float comp3 = 1 - u + 2 * (u * u);
    float comp4 = 1 + u - 2 * (u * u);
    float t = i * ((float)db / (float)da);

    // theoretical lower and upper bounds
    csc_int_t lb = (csc_int_t)floorf((comp3 / (comp1 * comp1 * comp2)) * t);
    csc_int_t ub = (csc_int_t)ceilf(((powf(comp2, 3) * comp4) / (comp1 * comp1)) * t);

    // starting point for the iteration
    csc_int_t sp = i - (ub - lb); // x coordinate
    csc_int_t y = lb; // y coordinate
    for (csc_int_t x = sp; x < i; x++) {
        // process overlined triangle down in Eq. (9) of the reference
        y += (2 * (x * db - y * da) >= 0);
    }
    *p_num_iter = ub - lb;
    return (D > A) ? (y + i) : y;
}

/**
 * \brief CSC based on the extended Bresenham's algorithm with practical bounds.
 *
 * \remarks For details, refer to the following paper:
 * - S. Kang and K. S. Kim, "Theoretical and practical bounds on the initial value
 *   of clock skew compensation algorithm immune to floating-point precision loss
 *   for resource-constrained wireless sensor nodes," Fiber and Integrated Optics,
 *   vol. 43, no. 3, pp. 111–121, Jun. 24, 2024.
 */
csc_int_t csc_eb3(const csc_int_t i, const csc_int_t D, const csc_int_t A, float epsilon, uint32_t *p_num_iter)
{
#ifndef CSC_NO_DIV_CHECK
    if (A == 0) {
        return 0;
    }
#endif
    csc_int_t da = A;
    csc_int_t db = (D > A) ? (D - A) : D; // Delta a and Delta b

    float u = powf(2, -24);
    float comp1 = 1 - u;
    float comp2 = 1 + 2 * u;
    float t = i * ((float)db / (float)da);

    // theoretical lower and upper bounds
    csc_int_t lb = (csc_int_t)floorf((comp1 / powf(comp2, 3))*t);
    csc_int_t ub = (csc_int_t)ceilf(powf(comp2, 3) * t);

    // starting point for the iteration
    csc_int_t sp = i - (ub - lb); // x coordinate
    csc_int_t y = lb; // y coordinate
    for (csc_int_t x = sp; x < i; x++) {
        // process overlined triangle down in Eq. (9) of the reference
        y += (2 * (x * db - y * da) >= 0);
    }
    *p_num_iter = ub - lb;
    return (D > A) ? (y + i) : y;
}

/**
 * \brief CSC based on the "direct search" algorithm.
 *
 * \remarks For details, refer to the following paper:
 * - K. S. Kim, "Direct search algorithm for clock skew compensation immune to
 *   floating-point precision loss," arXiv:2504.15039 [cs.NI], Apr. 2025.
 *   [Online]. Available: https://arxiv.org/abs/2504.15039
 */
csc_int_t csc_ds(const csc_int_t i, const csc_int_t D, const csc_int_t A, uint32_t *p_num_iter)
{
    csc_int_t j = 0;
    // csc_int_t k = floor(i*(float)D/(float)A + 0.5); // a starting point; floor() not working for uint64_t on TelosB platform
    csc_int_t k = (csc_int_t)(i*(float)D / (float)A + 0.5); // a starting point
    csc_int_t delta = (k - i) * A + i * (A - D); // "triangle down" to avoid overflow
    assert(delta == k * A - i * D); // for debugging

#ifdef CSC_NO_ITER_COUNT
    *p_num_iter = 1;
#else
    *p_num_iter = 0;
#endif
    if (delta == 0) {
        j = k;
#ifndef CSC_NO_ITER_COUNT
        (*p_num_iter)++;
#endif
    }
    else if (delta > 0) {
        while (true) {
            if (k == 0) {
                j = 0;
#ifndef CSC_NO_ITER_COUNT
                (*p_num_iter)++;
#endif
                break;
            }
            else {
                if (delta - A == 0) {
                    j = k - 1;
#ifndef CSC_NO_ITER_COUNT
                    (*p_num_iter)++;
#endif
                    break;
                }
                else if (delta - A > 0) {
                    k--;
                    delta -= A;
#ifndef CSC_NO_ITER_COUNT
                    (*p_num_iter)++;
#endif
                }
                else {
                    j = k - (ABS(delta - A) < ABS(delta)); // branchless programming
#ifndef CSC_NO_ITER_COUNT
                    (*p_num_iter)++;
#endif
                    break;
                }
            }
        } // end of while loop
    }
    else { // delta < 0
        while (true) {
            if (delta + A == 0) {
                j = k + 1;
#ifndef CSC_NO_ITER_COUNT
                (*p_num_iter)++;
#endif
                break;
            }
            else if (delta + A > 0) {
                j = k + (ABS(delta + A) < ABS(delta)); // branchless programming
#ifndef CSC_NO_ITER_COUNT
                (*p_num_iter)++;
#endif
                break;
            }
            else {
                k++;
                delta += A;
#ifndef CSC_NO_ITER_COUNT
                (*p_num_iter)++;
#endif
            }
        } // end of while loop
    } // delta < 0
    return j;
}

/**
 * \brief CSC based on the efficient implementation of the "direct search" algorithm
 *        without any loop and floating-point operation.
 *
 * \remarks For details, refer to the following paper:
 * - K. S. Kim, "Direct search algorithm for clock skew compensation immune to
 *   floating-point precision loss," arXiv:2504.15039 [cs.NI], Apr. 2025.
 *   [Online]. Available: https://arxiv.org/abs/2504.15039
 */
csc_int_t csc_ds2(const csc_int_t i, const csc_int_t D, const csc_int_t A)
{
    csc_int_t j = 0;
    //--------------------------------------------------------------------------
    // options for a starting point:
    // csc_int_t k = floor(i*(float)D/(float)A + 0.5); // a starting point; floor() not working for uint64_t on TelosB platform
    csc_int_t k = (csc_int_t)(i*(float)D/(float)A + 0.5); // a starting point
    // csc_int_t k = (i / A) * D; // a starting point without floating-point operations (but possible overflow when A > i)
    // csc_int_t k = i; // a starting point avoiding overflow resulting from the original one
    //--------------------------------------------------------------------------
    csc_int_t delta = (k - i) * A + i * (A - D); // "triangle down" to avoid overflow
    assert(delta == k * A - i * D); // for debugging

    if (delta == 0) {
        j = k;
    }
    else if (delta > 0) {
        k -= (delta / A);
        delta %= A;
        if (delta == 0) {
            j = k;
        }
        else {
            j = k - (ABS(delta - A) < ABS(delta)); // branchless programming
        }
    }
    else { // delta < 0
        k += (-delta / A);
        delta = delta % A; // N.B.: different from mathematical modulo
        if (delta + A > 0) {
            j = k + (ABS(delta + A) < ABS(delta)); // branchless programming
        }
        else {
            j = k;
        }
    } // delta < 0
    return j;
}

/**
 * \brief CSC based on the efficient implementation of the "direct search"
 * algorithm without any loop and floating-point operation but with an offset
 * from a prior run.
 *
 * \note This is an internal function and not to be directly called by users.
 *
 * \remarks For details, refer to the following paper:
 * - K. S. Kim, "Decomposition of integer division rounded to the nearest
 *   integer: Space-time trade-off in clock skew compensation," to be submitted
 *   to IEEE Trans. Singal Process., Apr. 2026."
 *   [Online]. Available: TBD
 */
csc_int_t csc_ds2_offset(const csc_int_t i, const csc_int_t D, const csc_int_t A, const csc_int_t offset, csc_int_t *remainder)
{
    short sign = 1;
    csc_int_t j = 0;
    //--------------------------------------------------------------------------
    // options for a starting point:
    // csc_int_t k = floor(i*(float)D/(float)A + 0.5); // a starting point; floor() not working for uint64_t on TelosB platform
    // csc_int_t k = (csc_int_t)(i*(float)D/(float)A + 0.5); // a starting point
    // csc_int_t k = (i / A) * D; // a starting point without floating-point operations (but possible overflow when A > i)
    csc_int_t k = i; // a starting point avoiding overflow resulting from the original one
    //--------------------------------------------------------------------------
    // "triangle down" with the offset to avoid overflow
#if DEBUG == 1
    csc_int_t delta;
    bool is_overflow;
    // is_overflow = INT_SUB_OVERFLOW_P(k, i) || INT_MUL_OVERFLOW_P(k - i, A) ||
    //     INT_SUB_OVERFLOW_P(A, D) || INT_MUL_OVERFLOW_P(i, A - D) ||
    //     INT_ADD_OVERFLOW_P((k - i) * A, i * (A - D)) ||
    //     INT_ADD_OVERFLOW_P((k - i) * A + i * (A - D), offset);
    is_overflow = INT_SUB_OVERFLOW_P(A, D) || INT_MUL_OVERFLOW_P(i, A - D) || INT_ADD_OVERFLOW_P(i * (A - D), offset);
    if (is_overflow) {
        fprintf(stderr, "Overflow in csc_ds2_offset: i=%"CSC_INT_PRI", D=%"CSC_INT_PRI", A=%"CSC_INT_PRI", offset=%"CSC_INT_PRI"\n",
            i, D, A, offset);
        // exit(EXIT_FAILURE);
        return -100; // return a negative number to indicate an overflow
    }
    else {
        // delta = (k - i) * A + i * (A - D) + offset;
        delta = i * (A - D) + offset;
    }
#else
    // csc_int_t delta = (k - i) * A + i * (A - D) + offset;
    csc_int_t delta = i * (A - D) + offset;
#endif
    assert(delta == k * A - i * D + offset); // for debugging

    if (delta == 0) {
        j = k;
        *remainder = delta;
    }
    else {
        k -= (delta / A);
        delta %= A;
        if (ABS(ABS(delta) - A) < ABS(delta)) {
            sign = (delta < 0) - (delta > 0); // sign of '-delta'; branchless programming
            j = k + sign;
            *remainder = delta + sign * A;
        }
        else {
            j = k;
            *remainder = delta;
        }
    }
    return j;
}

/**
 * \brief CSC based on the efficient implementation of the "direct search"
 * algorithm without any loop and floating-point operation but with an offset
 * from a prior run.
 *
 * \note This is an internal function and not to be directly called by users.
 *
 * \note Update: This version turns out to be slower than csc_ds2_offset due to
 * the additional arithmetic operations related with the "condition".
 *
 * \remarks For details, refer to the following paper:
 * - K. S. Kim, "Decomposition of integer division rounded to the nearest
 *   integer: Space-time trade-off in clock skew compensation," to be submitted
 *   to IEEE Trans. Singal Process., Apr. 2026."
 *   [Online]. Available: TBD
 */
csc_int_t csc_ds2_offset2(const csc_int_t i, const csc_int_t D, const csc_int_t A, const csc_int_t offset, csc_int_t *remainder)
{
    short condition, sign;
    // csc_int_t j = 0;
    //--------------------------------------------------------------------------
    // options for a starting point:
    // csc_int_t k = floor(i*(float)D/(float)A + 0.5); // a starting point; floor() not working for uint64_t on TelosB platform
    // csc_int_t k = (csc_int_t)(i*(float)D/(float)A + 0.5); // a starting point
    // csc_int_t k = (i / A) * D; // a starting point without floating-point operations (but possible overflow when A > i)
    csc_int_t k = i; // a starting point avoiding overflow resulting from the original one
    //--------------------------------------------------------------------------
    // "triangle down" with the offset to avoid overflow
#if DEBUG == 1
    csc_int_t delta;
    bool is_overflow;
    // is_overflow = INT_SUB_OVERFLOW_P(k, i) || INT_MUL_OVERFLOW_P(k - i, A) ||
    //     INT_SUB_OVERFLOW_P(A, D) || INT_MUL_OVERFLOW_P(i, A - D) ||
    //     INT_ADD_OVERFLOW_P((k - i) * A, i * (A - D)) ||
    //     INT_ADD_OVERFLOW_P((k - i) * A + i * (A - D), offset);
    is_overflow = INT_SUB_OVERFLOW_P(A, D) || INT_MUL_OVERFLOW_P(i, A - D) || INT_ADD_OVERFLOW_P(i * (A - D), offset);
    if (is_overflow) {
        fprintf(stderr, "Overflow in csc_ds2_offset: i=%"CSC_INT_PRI", D=%"CSC_INT_PRI", A=%"CSC_INT_PRI", offset=%"CSC_INT_PRI"\n",
            i, D, A, offset);
        // exit(EXIT_FAILURE);
        return -100; // return a negative number to indicate an overflow
    }
    else {
        // delta = (k - i) * A + i * (A - D) + offset;
        delta = i * (A - D) + offset;
    }
#else
    // csc_int_t delta = (k - i) * A + i * (A - D) + offset;
    csc_int_t delta = i * (A - D) + offset;
#endif
    assert(delta == k * A - i * D + offset); // for debugging

    if (delta == 0) {
        *remainder = delta;
        return k; // j = k
    }
    else {
        k -= (delta / A);
        delta %= A;
        condition = ABS(ABS(delta) - A) < ABS(delta); // condition for Cases 2.1 and 3.1
        sign = (delta < 0) - (delta > 0); // sign of '-delta'; branchless programming
        *remainder = delta + condition * sign * A;
        return k + condition * sign; // j = k + condition * sign
    }
}

/**
 * \brief CSC based on the additive decomposition of direct search (ADDS).
 *
 * \remarks For details, refer to the following paper:
 * - K. S. Kim, "Decomposition of integer division rounded to the nearest
 *   integer: Space-time trade-off in clock skew compensation," to be submitted
 *   to IEEE Trans. Singal Process., Apr. 2026."
 *   [Online]. Available: TBD
 */
csc_int_t csc_adds(const csc_int_t i, const csc_int_t D, const csc_int_t A, int N)
{
    csc_int_t i_n, i_n_last, j;
    csc_int_t offset, remainder = 0;

    i_n = i / N;
    i_n_last = i - (N - 1)*i_n;

    j = 0;
    offset = 0;
    for (int n = 0; n < N - 1; n++) {
        j += csc_ds2_offset(i_n, D, A, offset, &remainder);
        offset = remainder;
    }
    j += csc_ds2_offset(i_n_last, D, A, offset, &remainder);
    return j;
}

/**
 * \brief ADDS without offset for comparison.
 */
csc_int_t csc_adds_wo(const csc_int_t i, const csc_int_t D, const csc_int_t A, int N)
{
    csc_int_t i_n, i_n_last, j;
    csc_int_t offset, remainder;

    i_n = i / N;
    i_n_last = i - (N - 1)*i_n;

    j = 0;
    offset = 0;
    for (int n = 0; n < N - 1; n++) {
        j += csc_ds2_offset(i_n, D, A, offset, &remainder);
        // offset = remainder;
    }
    j += csc_ds2_offset(i_n_last, D, A, offset, &remainder);
    return j;
}

/**
 * \brief CSC based on the additive decomposition of direct search (ADDS).
 *
 * \note This function calls csc_ds2_offset2().
 *
 * \remarks For details, refer to the following paper:
 * - K. S. Kim, "Decomposition of integer division rounded to the nearest
 *   integer: Space-time trade-off in clock skew compensation," to be submitted
 *   to IEEE Trans. Singal Process., Apr. 2026."
 *   [Online]. Available: TBD
 */
csc_int_t csc_adds2(const csc_int_t i, const csc_int_t D, const csc_int_t A, int N)
{
    csc_int_t i_n, i_n_last, j;
    csc_int_t offset, remainder = 0;

    i_n = i / N;
    i_n_last = i - (N - 1)*i_n;

    j = 0;
    offset = 0;
    for (int n = 0; n < N - 1; n++) {
        j += csc_ds2_offset2(i_n, D, A, offset, &remainder);
        offset = remainder;
    }
    j += csc_ds2_offset2(i_n_last, D, A, offset, &remainder);
    return j;
}

/**
 * \brief A wrapper function for CSC based on the additive decomposition of direct search (ADDS).
 *
 * \remarks For details, refer to the following paper:
 * - K. S. Kim, "Decomposition of integer division rounded to the nearest
 *   integer: Space-time trade-off in clock skew compensation," to be submitted
 *   to IEEE Trans. Singal Process., Apr. 2026."
 *   [Online]. Available: TBD
 */
csc_int_t csc_adds_wrapper(const csc_int_t i, const csc_int_t D, const csc_int_t A)
{
    csc_int_t j;
    int N;
#if CSC_INT_SIZE == 4
    // if (i <= 10000000) { // 1E7
    if (i <= 100000) { // 1E5
        N = 1;
    }
    else {
        // N = i / 10000000;
        N = i / 100000;
    }
#elif CSC_INT_SIZE == 8
    // if (i <= 10000000000000) { // 1E13
    if (i <= 10000000000) { // 1E10
        N = 1;
    }
    else {
        // N = i / 10000000000000;
        N = i / 10000000000;
    }
#endif
    j = csc_adds(i, D, A, N);
    return j;
}

/**
 * \brief A wrapper function for CSC based on the additive decomposition of direct search (ADDS).
 *
 * \note This function calls csc_adds2().
 *
 * \remarks For details, refer to the following paper:
 * - K. S. Kim, "Decomposition of integer division rounded to the nearest
 *   integer: Space-time trade-off in clock skew compensation," to be submitted
 *   to IEEE Trans. Singal Process., Apr. 2026."
 *   [Online]. Available: TBD
 */
csc_int_t csc_adds2_wrapper(const csc_int_t i, const csc_int_t D, const csc_int_t A)
{
    csc_int_t j;
    int N;
#if CSC_INT_SIZE == 4
    // if (i <= 10000000) { // 1E7
    if (i <= 100000) { // 1E5
        N = 1;
    }
    else {
        // N = i / 10000000;
        N = i / 100000;
    }
#elif CSC_INT_SIZE == 8
    // if (i <= 10000000000000) { // 1E13
    if (i <= 10000000000) { // 1E10
        N = 1;
    }
    else {
        // N = i / 10000000000000;
        N = i / 10000000000;
    }
#endif
    j = csc_adds2(i, D, A, N);
    return j;
}

/**
 * \brief CSC based on the "Multiplicative Decomposition of Integer Division Rounded
 *        Division (MDID)" based on unsigned integers.
 */
csc_int_t csc_mdid(const csc_int_t i, const csc_int_t D, const csc_int_t A)
{
    csc_int_t a, b, quotient, remainder, rounding;
#ifndef CSC_NO_DIV_CHECK
    if (A == 0) {
        return 0;
    }
#endif
    if (i > D) {
        a = i;
        b = D;
    } else {
        a = D;
        b = i;
    }
#if DEBUG == 1
    bool is_overflow;
    is_overflow = INT_MUL_OVERFLOW_P(a / A, b);
    if (is_overflow) {
        fprintf(stderr, "Overflow in csc_mdid-quotient: i=%"CSC_INT_PRI", D=%"CSC_INT_PRI", A=%"CSC_INT_PRI"\n",
            i, D, A);
        // exit(EXIT_FAILURE);
        return -100; // return a negative number to indicate an overflow
    }
    else {
        quotient = (a / A) * b;
    }
#else
    quotient = (a / A) * b;
#endif
    remainder = a % A;
#if DEBUG == 1
    is_overflow = INT_MUL_OVERFLOW_P(remainder, b) || INT_ADD_OVERFLOW_P(remainder * b, A / 2);
    if (is_overflow) {
        fprintf(stderr, "Overflow in csc_mdid-rounding: i=%"CSC_INT_PRI", D=%"CSC_INT_PRI", A=%"CSC_INT_PRI"\n",
            i, D, A);
        // exit(EXIT_FAILURE);
        return -100; // return a negative number to indicate an overflow
    }
    else {
        rounding = (remainder * b + A / 2) / A;
    }
#else
    rounding = (remainder * b + A / 2) / A;
#endif
    return quotient + rounding;
}
