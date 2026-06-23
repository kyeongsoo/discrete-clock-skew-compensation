/**
 * \file csc-adds-ana.c
 * Investigation of the effects of the offset (i.e., \td_i) transfer
 * in the additive decomposition of direct search (ADDS) algorithm.
 *
 * \author Kyeong Soo (Joseph) Kim <kyeongsoo.kim@gmail.com>
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
    csc_int_t i = 1000000; // 1E6
    csc_int_t D = 1000000; // 1E6; corresponding to 1 s for 1 us clock resolution
    int Ns[] = {1, 10, 100, 1000};
    int u_max = 1000; // upper bound for a uniformly-distributed random integer in [1, u_max]
    char int_type_str[] = "int32";
#elif CSC_INT_SIZE == 8
    csc_int_t i = 1000000000000; // 1E12
    csc_int_t D = 1000000000000; // 1E12
    int Ns[] = {1, 10, 100, 1000};
    int u_max = 1000000; // upper bound for a uniformly-distributed random integer in [1, u_max]
    char int_type_str[] = "int64";
#endif
    int N_Ns = sizeof(Ns) / sizeof(Ns[0]);
    int N_samples = 1000; // number of samples for D
    int random_seed = 20260420; // for reproducibility
    csc_int_t A, b512_result, csc_abs_err, csc_err, j;
    csc_int_t csc_errs[N_samples] = {};
    long int u;

    // CSC algosrithms
    csc_int_t (*csc_algs[])(csc_int_t, csc_int_t, csc_int_t, int) = {
        csc_adds,
        csc_adds_wo
    };
    char *alg_names[] = {
        "Additive decomposition of direct search (ADDS)",
        "Additive decomposition of direct search (ADDS) w/o offset",
    };
    int N_algs = sizeof(csc_algs) / sizeof(csc_algs[0]);
    csc_int_t csc_err_min[N_algs], csc_err_max[N_algs];
    double csc_abs_err_sum[N_algs], csc_err_sum[N_algs];

    // output experimental settings
    printf("# Input parameters\n");
    printf("- i=%.0e\n", (float)i);
    printf("- D=%.0e\n", (float)D);
    printf("- u_max=%d\n", u_max);
    printf("- num_samples=%.0e\n", (float)N_samples);
    printf("- random_seed=%d\n", random_seed);
    printf("\n");
    printf("# Clock skew compensation results\n");

    srand48(random_seed); // seed for random number generator

    for (int n = 0; n < N_Ns; n++) {
        int N = Ns[n];

        // initialize statistics
        for (int a = 0; a < N_algs; a++) {
            csc_err_min[a] = CSC_INT_MAX;
            csc_err_max[a] = CSC_INT_MIN;
            csc_err_sum[a] = 0.0;
            csc_abs_err_sum[a] = 0.0;
        } // end of for() for a

        for (int k = 0; k < N_samples; k++) {
            u = (lrand48() % u_max) + 2;
            A = (csc_int_t) (D + u);
            b512_result = csc_b512(i, D, A);
            for (int a = 0; a < N_algs; a++) {
                j = csc_algs[a](i, D, A, N);

                // CSC error
                csc_err = b512_result - j;
                if (a == 1) {
                    csc_errs[k] = csc_err; // for preprocessing in Python
                }
                csc_err_sum[a] += csc_err;
                csc_abs_err = ABS(csc_err); // absolute error
#if DEBUG == 1
                if ((csc_abs_err > 1.0) && (a == 2)) {
                    printf("DEBUG: i=%.0e:\n", (float)i);
                    printf("DEBUG: D=%"CSC_INT_PRI"\n", D);
                    printf("DEBUG: A=%"CSC_INT_PRI"\n", A);
                    printf("DEBUG: j=%"CSC_INT_PRI"\n", j);
                    printf("DEBUG: b512_result=%"CSC_INT_PRI"\n", b512_result);
                    printf("DEBUG: CSC error=%"CSC_INT_PRI"\n", csc_err);
                    exit(1);
                }
#endif
                csc_abs_err_sum[a] += csc_abs_err;
                if (csc_err < csc_err_min[a]) {
                    csc_err_min[a] = csc_err;
                }
                if (csc_err > csc_err_max[a]) {
                    csc_err_max[a] = csc_err;
                }
            } // end of for() for a
        } // end of for() for k

        // store csc_errs for preprocessing in Python
        char filename[256];
        snprintf(filename, sizeof(filename), "./out/csc-adds-ana_%s_N%.0e.bin", int_type_str, (float)N);
        FILE *fp = fopen(filename, "wb");
        if (fp) {
            fwrite(csc_errs, sizeof(csc_err), N_samples, fp);
            fclose(fp);
        }
        else {
            fprintf(stderr, "Error opening file %s for writing\n", filename);
            exit(1);
        }

        // output CSC results
        printf("## N=%.0e:\n", (float)N);
        // printf("## N=%.0e:\n", (float)N);
        for (int a = 0; a < N_algs; a++) {
            printf("### %s\n", alg_names[a]);
            printf("- CSC error (min./max./avg.): %"CSC_INT_PRI"/%"CSC_INT_PRI"/%.4e\n",
                csc_err_min[a], csc_err_max[a], csc_err_sum[a]/N_samples);
            printf("- CSC mean absolute error: %.4e\n", csc_abs_err_sum[a]/N_samples);
            fflush(stdout);
        } // end of for() for a
    } // end of for() for N

    return 0;
}
