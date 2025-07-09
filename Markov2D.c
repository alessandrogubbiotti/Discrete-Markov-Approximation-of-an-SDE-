/* Markov2D_plot_ready.c */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define IDX(x, y) ((y) * Nx + (x))

int Nx =1000; 
int  Ny = 1000; 
int  N_particles = 100;  
int T_max = 10000000;
int snapshot_interval = 100000;

double beta = 1 ;

// Allocate flat arrays
double  *V, *mu, *gradVx, *gradVy, *cx, *cy, *mucx, *mucy, *divergence;


typedef struct{
	int time; 
	int particle_id; 
	int x; 
	int y; 
} Particle_State;

Particle_State *particles;

// GSL RNG setup
gsl_rng *rng;


// Configuration reader
//void read_config(const char *filename) {
//    FILE *fp = fopen(filename, "r");
//   if (!fp) { perror("Config file"); exit(1); }
//    fscanf(fp, "%d %d", &Nx, &Ny);
//    fscanf(fp, "%d", &N_particles);
//    fscanf(fp, "%d", &T_max);
//    fscanf(fp, "%lf", &beta);
//    fscanf(fp, "%d", &snapshot_interval);
//    fclose(fp);
//}

// User-defined potential V
//
void initialize_V(){
	for(int i = 0; i < Nx; i++){
		for( int j = 0; j < Ny ; j++){
		V[IDX(i,j)] = (abs(i - 10) *abs(j -10) + abs(i-50)*abs(j -60))/1000;
		}				
	}
}




//void initialize_V(double *V ){
//	for(int i = 0; i < Nx*Ny; i++) V[i] = 1; 
//
//    for (int i = 0; i < Nx * Ny; i++) {
//	}
//}

void compute_mu() {
    double Vmin = V[0];
    for (int i = 1; i < Nx * Ny; i++) {
        if (V[i] < Vmin) Vmin = V[i];
    }

    for (int i = 0; i < Nx * Ny; i++) {
        mu[i] = exp(-beta * (V[i] - Vmin));  // Numerically stable
    }

}

// Drift field (example: zero)
void drift(int x, int y, double *bx, double *by) {
    *bx = 0.0;
    *by = 0.0;
}


void normalize_mu() {
    double Z = 0.0;
    for (int i = 0; i < Nx * Ny; i++) Z += mu[i];
    if (Z == 0.0) {
        fprintf(stderr, "ERROR: Normalization factor Z is zero. Possibly overflow/underflow in mu.\n");
        exit(1);
    }
    for (int i = 0; i < Nx * Ny; i++) mu[i] /= Z;
}

void compute_fields() {
    compute_mu(V);
    for (int y = 0; y < Ny; y++) {
        for (int x = 0; x < Nx; x++) {
            int xp = (x + 1) % Nx, xm = (x - 1 + Nx) % Nx;
            int yp = (y + 1) % Ny, ym = (y - 1 + Ny) % Ny;

            double dx = 0.5 * (V[IDX(xp,y)] - V[IDX(xm,y)]);
            double dy = 0.5 * (V[IDX(x,yp)] - V[IDX(x,ym)]);
            gradVx[IDX(x,y)] = dx;
            gradVy[IDX(x,y)] = dy;

            double bx, by;
            drift(x, y, &bx, &by);
            cx[IDX(x,y)] = bx + dx;
            cy[IDX(x,y)] = by + dy;

            mu[IDX(x,y)] = exp(-beta * V[IDX(x,y)]);
        }
    }
    normalize_mu();
    for (int i = 0; i < Nx * Ny; i++) {
        mucx[i] = mu[i] * cx[i];
        mucy[i] = mu[i] * cy[i];
        divergence[i] = gradVx[i] * cx[i] + gradVy[i] * cy[i];
    }
}

void write_flat_csv(const char *filename, double *data) {
    FILE *fp = fopen(filename, "w");
    for (int y = 0; y < Ny; y++) {
        for (int x = 0; x < Nx; x++) {
            fprintf(fp, "%lf", data[IDX(x,y)]);
            if (x < Nx - 1) fprintf(fp, ",");
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}
//
//void write_positions_csv(const char *filename) {
//    FILE *fp = fopen(filename, "w");
//    for (int i = 0; i < N_particles; i++)
//        fprintf(fp, "%d,%d\n", x_particles[i], y_particles[i]);
//    fclose(fp);
//}
//
void simulate() {
    for (int i = 0; i < N_particles; i++) {
        particles[i].x = gsl_rng_uniform_int(rng, Nx);
        particles[i].y = gsl_rng_uniform_int(rng, Ny);
    }
FILE *traj_file = fopen("trajectory.bin", "wb");
	if(!traj_file) {
		perror(" Failed to opent the binary file\n");
		free(particles); 
		return; 
	}
	for (int t = 0; t <= T_max; t++) {
        for (int i = 0; i < N_particles; i++) {
		particles[i].time = t; 
            int dir = gsl_rng_uniform_int(rng, 4);
            int x = particles[i].x;
            int y = particles[i].y; 
            int xn = x, yn = y;
            if (dir == 0) xn = (x + 1) % Nx;
            else if (dir == 1) xn = (x - 1 + Nx) % Nx;
            else if (dir == 2) yn = (y + 1) % Ny;
            else yn = (y - 1 + Ny) % Ny;

            double V_old = V[IDX(x, y)];
            double V_new = V[IDX(xn, yn)];
            double dE = V_new - V_old;
            if (gsl_rng_uniform(rng) < exp(-beta * dE)) {
                particles[i].x = xn;
                particles[i].y  = yn;
            }
         fwrite( particles, sizeof(Particle_State), N_particles, traj_file); 
        }

     //   if (t % snapshot_interval == 0) {
       //     char fname[64];
         //   sprintf(fname, "positions_t%04d.csv", t);
           // write_positions_csv(fname);

            //int *hist = calloc(Nx * Ny, sizeof(int));
            //for (int i = 0; i < N_particles; i++) {
             //   hist[IDX(x_particles[i], y_particles[i])]++;
           // }
            //double *empirical = malloc(Nx * Ny * sizeof(double));
            //for (int i = 0; i < Nx * Ny; i++) empirical[i] = (double) hist[i] / N_particles;
           // sprintf(fname, "empirical_t%04d.csv", t);
           // write_flat_csv(fname, empirical);
           // free(hist); free(empirical);
        }
 

    fclose(traj_file); 
}

int main() {
 //   read_config("config.txt");
    gsl_rng_env_setup();
    rng = gsl_rng_alloc(gsl_rng_mt19937);
    V = malloc(Nx * Ny * sizeof(double));
    mu = malloc(Nx * Ny * sizeof(double));
    initialize_V();
    compute_mu(); 
    gradVx = malloc(Nx * Ny * sizeof(double));
    gradVy = malloc(Nx * Ny * sizeof(double));
    cx = malloc(Nx * Ny * sizeof(double));
    cy = malloc(Nx * Ny * sizeof(double));
    mucx = malloc(Nx * Ny * sizeof(double));
    mucy = malloc(Nx * Ny * sizeof(double));
    divergence = malloc(Nx * Ny * sizeof(double));
    particles = malloc(N_particles * sizeof(Particle_State));
        
    compute_fields();
    write_flat_csv("V.csv", V);
    write_flat_csv("mu.csv", mu);
    write_flat_csv("gradVx.csv", gradVx);
    write_flat_csv("gradVy.csv", gradVy);
    write_flat_csv("c_x.csv", cx);
    write_flat_csv("c_y.csv", cy);
    write_flat_csv("muc_x.csv", mucx);
    write_flat_csv("muc_y.csv", mucy);
    write_flat_csv("divV_dot_c.csv", divergence);
    simulate();
    gsl_rng_free(rng);
    free(V); free(mu); free(gradVx); free(gradVy);
    free(cx); free(cy); free(mucx); free(mucy); free(divergence);
    free(particles);
    return 0;
}

