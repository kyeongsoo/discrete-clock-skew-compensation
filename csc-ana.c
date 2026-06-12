/**
 * \file csc-ana.c
 * Analysis of clock skew compensation (CSC) algorithms
 * based on integer linear scaling rounded to the nearest integer
 *
 * \author Kyeong Soo (Joseph) Kim <kyeongsoo.kim@gmail.com>
 *
 * \remarks The results are published in the following paper:
 * - Kyeong Soo Kim, "Space-time trade-off in integer linear scaling
 *   rounded to the nearest integer through multiplicative and additive
 *   decomposition," arXiv e-prints arXiv:2605.21400v [cs.DS], May 2026.
 *   [Online]. Available: https://arxiv.org/abs/2605.21400
 */

#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "csc.h"
#include "drand48.h"

int main(int argc, char **argv)
{
    // experimental parameters
#if CSC_INT_SIZE == 4
    // csc_int_t D = 1000000; // 1E6; corresponding to 1 s for 1 us clock resolution
    csc_int_t D = 100000000; // 1E8; corresponding to 100 s for 1 us clock resolution
    csc_int_t is[] = {
        // 1000, // 1E3
        // 10000,
        100000, // 1E5
        1000000, // 1E6
        10000000,
        100000000,
        1000000000 // 1E9
    };
#elif CSC_INT_SIZE == 8
    // csc_int_t D = 1000000000; // 1E9; corresponding to 1 s for 1 ns clock resolution
    csc_int_t D = 1000000000000; // 1E12; corresponding to 1E3 s for 1 ns clock resolution
    csc_int_t is[] = {
        // 1000000, // 1E6
        // 10000000,
        // 100000000,
        // 1000000000, // 1E9
        10000000000, // 1E10
        100000000000,
        1000000000000, // 1E12
        10000000000000,
        100000000000000,
        1000000000000000,
        10000000000000000,
        100000000000000000,
        1000000000000000000 // 1E18
    };
#elif CSC_INT_SIZE == 16
    // TODO: Set D and is[] for 128-bit integers
    csc_int_t D = 1000000000000000000000; // 1E21; corresponding to 1E3 s for 1 as (atto second=1E-18 s) clock resolution
    csc_int_t is[] = {
        1000000 // 1E38
    };
#endif
    int N_is = sizeof(is) / sizeof(is[0]);
    // int N_samples = 1000000; // number of samples for D
    int N_samples = 100; // number of samples for D
    int random_seed = 20260421; // for reproducibility
    int skew_max = 100; // skew bound in ppm
    csc_int_t A, b512_result, csc_abs_err, csc_err, i, j;
    long int skew;


    // CSC algosrithms
    csc_int_t (*csc_algs[])(csc_int_t, csc_int_t, csc_int_t) = {
        csc_b256,
        csc_b128,
        csc_dp,
        csc_sp,
        csc_mdid,
        csc_adds_wrapper
    };
    char *alg_names[] = {
        "Binary256 FP division",
        "Binary128 FP division",
        "Double-precision FP division",
        "Single-precision FP division",
        "Multiplicative decomposition of integer division (MDID)",
        "Additive decomposition of direct search (ADDS)"
    };
    int N_algs = sizeof(csc_algs) / sizeof(csc_algs[0]);
    csc_int_t csc_err_min[N_algs], csc_err_max[N_algs];
    double csc_abs_err_sum[N_algs], csc_err_sum[N_algs];

    // output experimental settings
    printf("# Input parameters\n");
    printf("- D=%.0e\n", (float)D);
    printf("- skew_max=%d\n", skew_max);
    printf("- num_samples=%.0e\n", (float)N_samples);
    printf("- random_seed=%d\n", random_seed);
    printf("\n");
    printf("# Clock skew compensation results\n");

    srand48(random_seed); // seed for random number generator

    for (int n = 0; n < N_is; n++) {
        i = is[n];

        // initialize statistics
        for (int a = 0; a < N_algs; a++) {
            csc_err_min[a] = CSC_INT_MAX;
            csc_err_max[a] = CSC_INT_MIN;
            csc_err_sum[a] = 0.0;
            csc_abs_err_sum[a] = 0.0;
        } // end of for() for a

        for (int a = 0; a < N_algs; a++) {
            for (int k = 0; k < N_samples; k++) {
                skew = (lrand48() % (2*skew_max + 1)) - skew_max; // in ppm
                A = (csc_int_t) ((1 + skew*1.0E-6)*D);
                b512_result = csc_b512(i, D, A);
                j = csc_algs[a](i, D, A);

                // handle overflow
                if (j == -100) {
                    fprintf(stderr, "Overflow in %s: i=%"CSC_INT_PRI", D=%"CSC_INT_PRI", A=%"CSC_INT_PRI"\n",
                        alg_names[a], i, D, A);
                    csc_err_sum[a] = -100.0;
                    break; // skip this algorithm
                }

                // CSC error
                csc_err = b512_result - j;
                csc_err_sum[a] += csc_err;
                csc_abs_err = ABS(csc_err); // absolute error
#if DEBUG == 1
                if ((csc_abs_err > 1.0) && (a == 4)) {
                    printf("DEBUG: i=%.0e:\n", (float)i);
                    printf("DEBUG: D=%"CSC_INT_PRI"\n", D);
                    printf("DEBUG: A=%"CSC_INT_PRI"\n", A);
                    printf("DEBUG: j=%"CSC_INT_PRI"\n", j);
                    printf("DEBUG: b512_result=%"CSC_INT_PRI"\n", b512_result);
                    printf("DEBUG: CSC error=%"CSC_INT_PRI"\n", csc_err);
                    exit(1);
                }
#endif
                csc_abs_err_sum[a] += csc_abs_err;                if (csc_err < csc_err_min[a]) {
                    csc_err_min[a] = csc_err;
                }
                if (csc_err > csc_err_max[a]) {
                    csc_err_max[a] = csc_err;
                }
            } // end of for() for k (i.e., i samples)
        } // end of for() for a (i.e., algorithms)

        // output CSC results
        printf("## i=%.0e:\n", (float)i);
        for (int a = 0; a < N_algs; a++) {
            printf("### %s\n", alg_names[a]);
            if (csc_err_sum[a] == -100.0) {
                printf("- CSC error: Overflow occurred; results are not available.\n");
            }
            else {
                printf("- CSC error (min./max./avg.): %"CSC_INT_PRI"/%"CSC_INT_PRI"/%.4e\n",
                    csc_err_min[a], csc_err_max[a], csc_err_sum[a]/N_samples);
                printf("- CSC mean absolute error: %.4e\n", csc_abs_err_sum[a]/N_samples);
            }
            fflush(stdout);
        } // end of for() for a
    } // end of for() for i

    return 0;
}
