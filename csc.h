/**
 * \file csc.h
 * Analysis of clock skew compensation (CSC) algorithms.
 *
 * \author Kyeong Soo (Joseph) Kim <kyeongsoo.kim@gmail.com>
 *
 * \note The following options are controlled by macro definitions:
 * - CSC_INT_SIZE: The number of bytes for 'i', 'D', and 'A' (4 or 8).
 *
 * \remarks The results are published in the following paper:
 * - Kyeong Soo Kim, "Space-time trade-off in integer linear scaling
 *   rounded to the nearest integer through multiplicative and additive
 *   decomposition," arXiv e-prints arXiv:2605.21400v [cs.DS], May 2026.
 *   [Online]. Available: https://arxiv.org/abs/2605.21400
 *
 * \remarks Copyright (c) 2026 Kyeong Soo (Joseph) Kim
 *
 * \remarks SPDX-License-Identifier: MIT
 */

#ifndef CSC_H_
#define CSC_H_

// check overflow in integer operations (for debugging)
#if DEBUG == 1
#define INT_ADD_OVERFLOW_P(a, b) \
   __builtin_add_overflow_p (a, b, (__typeof__ ((a) + (b))) 0)
#define INT_SUB_OVERFLOW_P(a, b) \
   __builtin_sub_overflow_p (a, b, (__typeof__ ((a) - (b))) 0)
#define INT_MUL_OVERFLOW_P(a, b) \
   __builtin_mul_overflow_p (a, b, (__typeof__ ((a) * (b))) 0)
#endif

#if CSC_INT_SIZE == 4
typedef uint32_t csc_uint_t;
#define CSC_UINT_MAX UINT32_MAX
#define CSC_UINT_MIN 0
#define CSC_UINT_PRI PRIu32
typedef int32_t csc_int_t;
#define CSC_INT_MAX INT32_MAX
#define CSC_INT_MIN INT32_MIN
#define CSC_INT_PRI PRId32
#elif CSC_INT_SIZE == 8
typedef uint64_t csc_uint_t;
#define CSC_UINT_MAX UINT64_MAX
#define CSC_UINT_MIN 0
#define CSC_UINT_PRI PRIu64
typedef int64_t csc_int_t;
#define CSC_INT_MAX INT64_MAX
#define CSC_INT_MIN INT64_MIN
#define CSC_INT_PRI PRId64
#elif CSC_INT_SIZE == 16 && defined(__SIZEOF_INT128__)
typedef __uint128_t csc_uint_t;
#define CSC_UINT_MAX (~(__uint128_t)0)
#define CSC_UINT_MIN 0
#define CSC_UINT_PRI PRIu128
typedef __int128_t csc_int_t;
#define CSC_INT_MAX ((__int128_t)(( (__uint128_t)1 << 127 ) - 1))
#define CSC_INT_MIN (-CSC_INT_MAX - 1)
#define CSC_INT_PRI PRId128
#else
#error Unsupported CSC_INT_SIZE
#endif

// type-independent implementation
#ifndef ABS
#define ABS(x) (((x) < 0) ? -(x) : (x))
#endif

#if defined(__x86_64__) || defined(_M_X64) || defined(__amd64__)
csc_int_t csc_b128(const csc_int_t i, const csc_int_t D, const csc_int_t A);
csc_int_t csc_b256(const csc_int_t i, const csc_int_t D, const csc_int_t A);
csc_int_t csc_b512(const csc_int_t i, const csc_int_t D, const csc_int_t A);
#endif
csc_int_t csc_dp(const csc_int_t i, const csc_int_t D, const csc_int_t A);
csc_int_t csc_sp(const csc_int_t i, const csc_int_t D, const csc_int_t A);
csc_int_t csc_eb1(const csc_int_t i, const csc_int_t D, const csc_int_t A, float epsilon, uint32_t *p_num_iter);
csc_int_t csc_eb2(const csc_int_t i, const csc_int_t D, const csc_int_t A, float epsilon, uint32_t *p_num_iter);
csc_int_t csc_eb3(const csc_int_t i, const csc_int_t D, const csc_int_t A, float epsilon, uint32_t *p_num_iter);
csc_int_t csc_ds(const csc_int_t i, const csc_int_t D, const csc_int_t A, uint32_t *p_num_iter);
csc_int_t csc_ds2(const csc_int_t i, const csc_int_t D, const csc_int_t A);
csc_int_t csc_ds2_offset(const csc_int_t i, const csc_int_t D, const csc_int_t A, const csc_int_t offset, csc_int_t *remainder);
csc_int_t csc_ds2_offset2(const csc_int_t i, const csc_int_t D, const csc_int_t A, const csc_int_t offset, csc_int_t *remainder);
csc_int_t csc_adds(const csc_int_t i, const csc_int_t D, const csc_int_t A, int N);
csc_int_t csc_adds2(const csc_int_t i, const csc_int_t D, const csc_int_t A, int N);
csc_int_t csc_adds_wo(const csc_int_t i, const csc_int_t D, const csc_int_t A, int N);
csc_int_t csc_adds_wrapper(const csc_int_t i, const csc_int_t D, const csc_int_t A);
csc_int_t csc_adds2_wrapper(const csc_int_t i, const csc_int_t D, const csc_int_t A);
csc_int_t csc_mdid(const csc_int_t i, const csc_int_t D, const csc_int_t A);

#endif /* CSC_H_ */
