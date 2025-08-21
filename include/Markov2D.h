#ifndef MARKOV2D_H
#define MARKOV2D_H

#include <math.h>   // if you need fabs, etc.

/* ==== Structs ==== */
typedef struct {
    double north, south, east, west;
} Vector;

typedef struct {
    long double x;
    long double y;
} Form;

typedef struct {
    long double x;
    long double y;
} Activity;


// Particle_State gives the state of a particle (position)
typedef struct{
	int x; 
	int y; 
} Particle_State;

/* ==== Inline helper ==== */
static inline int IDX(int x, int y, int N) {
    x %= N;
    y %= N;
    if (x < 0) x += N;
    if (y < 0) y += N;
    return x * N + y;
}

inline long double squarel(long double v) {
    return v * v;
}

// Kahan compensated addition for long double
static inline void kahan_add(long double term, long double *sum, long double *comp) {
    long double y = term - *comp;
    long double t = *sum + y;
    *comp = (t - *sum) - y;
    *sum = t;
}

// Numerically-stable log-sum-exp *style* sum for exp(-beta V):
// shift by Vmin to avoid overflow/underflow; uses Kahan for the sum.
static long double stable_partition_sum(const double *V, double beta, int NN, double *Vmin_out) {
    double Vmin = V[0];
    for (int i = 1; i < NN; ++i) if (V[i] < Vmin) Vmin = V[i];
    *Vmin_out = Vmin;

    long double S = 0.0L, c = 0.0L;
    for (int i = 0; i < NN; ++i) {
        long double term = expl(-(long double)beta * ((long double)V[i] - (long double)Vmin));
        kahan_add(term, &S, &c);
    }
    return S; // Note: true Z = e^{-beta Vmin} * S, but we never need Z explicitly
}

#endif // MARKOV2D_H

