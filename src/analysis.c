#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "analysis.h"
#include "Markov2D.h"   // if you want access to IDX(), Form, etc.

// 1. Measure distances
void measure_distances(double *empirical, long double *m_beta, int N) {
    long double L1 = 0.0L, L2 = 0.0L, Linf = 0.0L;
    for (int i = 0; i < N*N; i++) {
        long double diff = fabsl((long double)empirical[i] - m_beta[i]);
        L1 += diff;
        L2 += diff * diff;
        if (diff > Linf) Linf = diff;
    }
    L2 = sqrtl(L2);

    printf("Distances between empirical and invariant measure:\n");
    printf("  L1 = %.12Lg\n", L1);
    printf("  L2 = %.12Lg\n", L2);
    printf("  Linf = %.12Lg\n", Linf);
}

// 2. Save heatmap differences
void save_measure_difference(const char *dirname, double *empirical, long double *m_beta, int N) {
    char filepath[512];
    snprintf(filepath, sizeof(filepath), "%s/measure_diff.txt", dirname);
    FILE *fp = fopen(filepath, "w");
    if (!fp) { perror("write diff"); return; }

    for (int x = 0; x < N; x++) {
        for (int y = 0; y < N; y++) {
            int idx = IDX(x,y,N);
            long double diff = fabsl((long double)empirical[idx] - m_beta[idx]);
            fprintf(fp, "%d %d %.8Lf\n", x, y, diff);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}

// 3. Estimate c from equilibrium current
void estimate_c(Form *c_est, Form *eq_current, long double *m_beta, double *V, double beta, int N) {
    for (int x = 0; x < N; x++) {
        for (int y = 0; y < N; y++) {
            int idx = IDX(x,y,N);
            int idx_east  = IDX(x, (y+1)%N, N);
            int idx_north = IDX((x+1)%N, y, N);

            long double denom_x = m_beta[idx] * expl(-0.5L * beta * (V[idx_east] - V[idx]));
            c_est[idx].x = (fabsl(denom_x) > 1e-14L) ? eq_current[idx].x / denom_x : 0.0L;

            long double denom_y = m_beta[idx] * expl(-0.5L * beta * (V[idx_north] - V[idx]));
            c_est[idx].y = (fabsl(denom_y) > 1e-14L) ? eq_current[idx].y / denom_y : 0.0L;
        }
    }
}


