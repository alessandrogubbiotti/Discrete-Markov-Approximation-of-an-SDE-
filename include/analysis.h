#ifndef ANALYSIS_H
#define ANALYSIS_H

#include "Markov2D.h"   // for Form, Vector, IDX, etc.

/**
 * Compute distances between the empirical measure and the invariant measure.
 * Distances are reported in L1, L2, and Linf norms.
 */
void measure_distances(double *empirical, long double *m_beta, int N);

/**
 * Save the local absolute difference between empirical and invariant measure
 * into a text file (for plotting as a heatmap).
 *
 * Output format per line: x y diff
 */
void save_measure_difference(const char *dirname,
                             double *empirical,
                             long double *m_beta,
                             int N);

/**
 * Estimate the antisymmetric field c from the equilibrium current.
 *
 * c_est[x,y].x corresponds to the east edge, c_est[x,y].y to the north edge.
 */
void estimate_c(Form *c_est,
                Form *eq_current,
                long double *m_beta,
                double *V,
                double beta,
                int N);

#endif // ANALYSIS_H

