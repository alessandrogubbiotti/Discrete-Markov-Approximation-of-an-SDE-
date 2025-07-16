/* Markov2D_plot_ready.c */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>

#define IDX(x, y) ((y) * Nx + (x))
int N = 101; 
int Nx =101; 
int  Ny = 101; 
int  N_particles =10000 ;  
int T_max = 1000000;
int snapshot_interval = 10000;

double beta = 0.5;

// Allocate flat arrays
double  *V, *c_x, *c_y, *mu, *gradVx, *gradVy, *cx, *cy, *mucx, *mucy, *divergence;
double *empirical_flat , *north_flat, *south_flat, *east_flat, *west_flat, *b_north_flat, *b_east_flat, *b_west_flat,  *b_south_flat;  

// Two dimensional versions of the above arrays



     //Stores the matrix versions of the above arrays. The flat versions are there for accessibility reasons in C
   	// INSTANTANEOUS EMPIRICAL  MEASURE 
	 double **empirical;
	// INSTANTANEOUS EMPIRICAL FLOW  
    double **east; 
    double **north;
    double **west; 
    double **south; 
    	// ESTIMAITON OF THE VECTOR FIELD
    double **b_north; 
    double **b_south;
    double **b_east; 
    double **b_west; 
      

typedef struct{
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



// The distance (discontinuous at the boundary )
double distance( int x, int y, int a, int b ){
	return ((x  -a)*(x- a  )+(y- b)*(y-b))/(Nx); 
}

// One Minimum
void initialize_V(){
	for(int i = 0; i < Nx; i++){
		for( int j = 0; j < Ny ; j++){
		V[IDX(i,j)] = distance(i,j , 30, 30 ) + distance(i, j,  70, 70 );
		}				
	}
}

// Distance (continuous at the boundary )
double distance2 (int x, int y, int a, int b ){
    	int dx = abs( x - a ); 
	int dy = abs(y - b); 
	dx = dx > Nx / 2 ? Nx - dx : dx;
   	dy = dy > Ny / 2 ? Ny - dy : dy;
	return dx + dy; 
}

void initialize2_V(){
	for(int i = 0; i < Nx; i++){
		for( int j = 0; j < Ny ; j++){
		V[IDX(i,j)] = (distance(i, j, 30, 30 ))/100;
		}				
	}
}
void initialize3_V(){
	for(int i = 0; i < Nx; i++){
		for( int j = 0; j < Ny ; j++){
		V[IDX(i,j)] = (abs(i - 10) *abs(j -10) + abs(i-50)*abs(j -60))/100;
		}				
	}
}

void initialize4_V() {

	for(int i = 0; i < Nx; i++){
		for( int j = 0; j < Ny ; j++){
		V[IDX(i,j)] = fmin( distance2(i,j, 30, 30 ) , distance2(i, j, 70 , 70))/100;
		}				
	}
}

//We take the rates (c^s(x,y) + c^a(x,y) ) e^{-\beta/2 V(y) - V(x)}, fo
// for an antisymmetric edge functiokn c^a- Usually c^s = 1, and c^s is small. 
// In order for e^{-\betaV}  to be the invariant measure, the following has to be satisfied 
// sum_{y \, y ~ x} c^a(x, y ) (e^{-\beta/2 V(y) } + e^{-\beta/2 V(x)}) = 0 
// The next function computes the above quantity 
// In the case of "orthogonal forces" this is clear, sicne c has itself divergence 0 and when $c$ is different from zero e^{-\beta/2  V(y)} + e^{-\beta/2 V} is constant.
void check_orthogonality() {
	double max = 0 ;
	double edge_east, edge_west, edge_north, edge_south, div;
	int i1, j1; 
	for(int i = 0; i < N; i++ ){
		for(int j = 0; j <N; j++){ 
			j1 = (j + 1)  % Nx; 
			edge_east =  c_x[ IDX(i,j)] *exp(-beta* V[IDX(i, j1)]  ); // Notice that the east edge is obtained by changing j 	
			j1 = (j - 1 + Nx) % Nx; 
			edge_west =  - c_x[ IDX(i,j1)] *exp(-beta* V[IDX(i, j1)]  );	

			i1 = (i + 1) % Ny;
			edge_north =  +  c_y[ IDX(i,j)] *exp(-beta* V[IDX(i1, j)]  );
			i1 = (i - 1 + Ny) % Ny;  
			edge_south =  - c_y[ IDX(i1,j)] *exp(-beta* V[IDX(i1, j)]  );	
			div = edge_east + edge_north + edge_west + edge_south;	
			
			if (div > max ) max = div; 
		}
	}
	printf(" The orthogonality check gives: %f\n", max); 
		
}


// Here we mimich the case in which $V$ and $c$ are radial and orthogonal, and $c$ is of divergence 0
void initialize5_V(){
int max = (N +1) / 2;
if (N % 2 != 1) {
        fprintf(stderr, "Error: N must be odd, but N = %d\n", N);
        return;  // exit with error code
    }	
// Define the N/2-1 different values of $V$
double *a = malloc(((N+1)/2)* sizeof(double) ); 
// Define the N/2 - 1 possible values fo $c$
double *b = malloc((N+1)/2  * sizeof(double)); 
	for(int i = 0; i < max; i++){
	a[i ] = (double)(i*i)/N; // It took a lot of time to understand that the casting went wrong  
	b[i] = sin (i); 
	}
	for(int i = 0; i < Nx* Ny; i++){
	c_x[i] = 0; 
	c_y[i ] = 0; 
	}
	V[IDX((N-1)/2 , (N - 1)/2) ] = a[0]; 
	for(int i = 1; i <  max; i++ ){
		int max2 = (N -1 )/2 +i ; 		
		for(int j = (N-1)/2 - i  ; j < max2; j++){

			V[IDX( (N-1)/2 + i , j )] = a[i];//Right side of the square  minus hte hifhrst point 
			V[IDX((N-1)/2 - i, N - j- 1 )  ] = a[i];  	// left side of the square minus the lowest point 
			V[IDX( N-1 - j, (N-1)/2 + i )] = a[i]; // upper sidei. The index that is boving is j 
			V[ IDX(  j, (N-1)/2 - i) ] = a[i]; // lower side 
		// c_x denotes the magnitude of the vector in the east dirsction, that is c((a,b ) , (a, b) + e_), while c_y in the north 
			c_x[IDX( (N-1)/2 + i , j )] = -b[i]; //upper  side of the square 
			c_x[IDX((N-1)/2 - i,j)  ] = b[i];  	// lower  side of the square
			c_y[IDX( j , (N-1)/2 + i )] =  b[i]; // irght side
			c_y[ IDX( j, (N-1)/2 - i) ] = -b[i];  // left side 
			

			}
		
	}
//			for (int i= N-1; i > -1; i--){
//				for (int j = 0; j <N ; j++ ){
//				printf("  (%f, %f)  " , c_x[IDX(i, j )], c_y[IDX(i,j)]); 
//				}
//				printf("\n"); 
//			}
//			for (int i= N-1; i > -1; i--){
//				for (int j =0; j < N; j++ ){
//				printf(" %f ", V[IDX(i,j)]); 
//				}
//				printf("\n"); 
//			}
}


void double_well(){
if (N % 4 != 1) {
        fprintf(stderr, "Error: N must be a multiple of 2 such which is not a multiple of 4, but N = %d\n", N);
        return;  // exit with error code
    }	
// Define the N/2-1 different values of $V$
double *a = malloc(((N-1)/2 + 1)* sizeof(double) ); 
// Define the N/2 - 1 possible values fo $c$
double *b = malloc(((N- 1)/2 + 1)  * sizeof(double));
double *b2 = malloc (((N-1 )/2+1) *sizeof(double));  
	for(int i = 0; i <= (N-1)/2 ; i++){
	a[i ] = (double)(i*i)/N+10; // It took a lot of time to understand that the casting went wrong  
	b[i] = i; 
	b2[i] = - i; 
	}
	for(int i = 0; i < Nx* Ny; i++){
	c_x[i] = 0; 
	c_y[i ] = 0; 
	}
	// Here we define hte potential at the two minima
	V[IDX((N-1)/2 , (N-1)/4 )] =  a[0]; 
	V[IDX((N - 1)/2 , 3*(N - 1)/4)] = a[0]; 
	//Here we loop onn the concentroc squares arond (N/4 -1, N/2 - 1)
	// The number of loops is i = 1, ..., N/4 - 1
	for(int i = 1; i <=   (N-1)/4  ; i++ ){ //Number of concentric squares 		
		for(int j = - i   ; j < i ; j++){ //side of the i-th concentric square
			
			V[IDX( (N-1)/2 + i  , (N-1)/4  -  j)] = a[i];//UPper side of the square 
			V[IDX( (N-1)/2  - i  , (N-1)/4 + j ) ] = a[i];//Lower side of the square  
			V[IDX( (N-1)/2   + j  , (N-1)/4 + i ) ] = a[i];//Right side of the square  min 
			V[IDX( (N-1)/2   -j  , (N-1)/4 - i ) ] = a[i];//Left side of the squaret 
		// c_x denotes the magnitude of the vector in the east dirsction, that is c((a,b ) , (a, b) + e_), while c_y in the north 
			c_x[IDX( (N-1)/2 +i  , (N-1)/4  + j  )] = -b[i]; //upper  side of the square 
			c_x[IDX((N-1)/2  - i  ,(N-1)/4  +j)  ] = b[i];  	// lower  side of the square
			c_y[IDX( (N-1)/2  +j     , (N-1)/4 + i  )] =  b[i]; // right side
			c_y[ IDX( (N-1)/2  +   j  , (N-1)/4 - i) ] = -b[i];  // left side 
			}
		
	}

	//Here we loop onn the vertical and horizontal lines,
	// The number of loops is i = 1, ..., N/4 - 1
	for(int i = 1; i <=  (N-1)/4  ; i++ ){ //Number of concentric squares 		
		for(int j = - i   ; j < i ; j++){ //side of the i-th concentric square
			printf("For the pair of indexes (%d,%d): \n"
       " Upper (%d, %d), Lower  (%d, %d),  Right (%d, %d),  Left(%d, %d)\n",
       i, j,
       (N-1)/2 + i, 3*(N-1)/4 - j,
       (N-1)/2 - i, 3*(N-1)/4 + j,
       (N-1)/2 +j , (N-1)/2 +i,
       3*(N-1)/4  + j, (N-1)/2 - i);
			V[IDX( (N-1)/2  +i , 3*(N-1)/4  -   j)] = a[i];//a[i];//upper side of the square 
			V[IDX((N-1)/2  - i  , 3*(N-1)/4   + j ) ] = a[i];//lower side of the square  
			V[IDX( (N-1)/2 + j  , 3*(N-1)/4 + i ) ] = a[i];//a[i];//right side of the square  min 
			V[IDX( (N-1)/2 - j  , 3*(N-1)/4 - i ) ] = a[i];//left side of the squaret 
		// c_x denotes the magnitude of the vector in the east dirsction, that is c((a,b ) , (a, b) + e_), while c_y in the north 
			c_x[IDX((N-1)/2+ i  , 3*(N-1)/4   +  j  )] = -b2[i]; //upper  side of the square 
			c_x[IDX((N-1)/2 - i  ,3*(N-1)/4  +j)  ] = b2[i];  	// lower  side of the square
			c_y[IDX( (N-1)/2  + j   , 3*(N-1)/4 + i )] =  b2[i]; // irght side
			c_y[ IDX( (N-1)/2 +   j , 3*(N-1)/4 - i) ] = -b2[i];  // left side 
			}
		
	}
	//Here we loop over  the external lines
	for(int i = 0; i <  (N-1)/4  ; i++ ){ // The  lines
		for(int j = 0   ; j < N ; j++){ // The vertexes in the line 
			
			V[IDX(  i, j    )] = a[ (N-1)/ 2   - i  ]; // Lower line  
			V[IDX(  N-1 -  i, j  ) ] =  a[(N -1)/2 - i ];  	 // upper line 
			c_x[IDX(  i, j  )] = b[ (N-1)/2 - i]; // The line at height i 
			c_x[IDX(  N-1 -  i , j)  ] = - b[(N-1)/2 - i];  	// The line at height (N-1 )/2 + (N-1)/4 + 1 + i 
		}
		
	}

		for (int i= N-1; i > -1; i--){
				for (int j = 0; j <N ; j++ ){
				printf("  (%f, %f)  " , c_x[IDX(i, j )], c_y[IDX(i,j)]); 
				}
				printf("\n"); 
			}
			for (int i= N-1; i > -1; i--){
				for (int j =0; j < N; j++ ){
				printf(" %f ", V[IDX(i,j)]); 
				}
				printf("\n"); 
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

//void compute_fields() {
//    compute_mu();
//    for (int y = 0; y < Ny; y++) {
//        for (int x = 0; x < Nx; x++) {
//            int xp = (x + 1) % Nx, xm = (x - 1 + Nx) % Nx;
//            int yp = (y + 1) % Ny, ym = (y - 1 + Ny) % Ny;
//
//            double dx = 0.5 * (V[IDX(xp,y)] - V[IDX(xm,y)]);
//            double dy = 0.5 * (V[IDX(x,yp)] - V[IDX(x,ym)]);
//            gradVx[IDX(x,y)] = dx;
//            gradVy[IDX(x,y)] = dy;
//
//            double bx, by;
//            drift(x, y, &bx, &by);
//            cx[IDX(x,y)] = bx + dx;
//            cy[IDX(x,y)] = by + dy;
//
//            mu[IDX(x,y)] = exp(-beta * V[IDX(x,y)]);
//        }
//    }
//    normalize_mu();
//    for (int i = 0; i < Nx * Ny; i++) {
//        mucx[i] = mu[i] * cx[i];
//        mucy[i] = mu[i] * cy[i];
//        divergence[i] = gradVx[i] * cx[i] + gradVy[i] * cy[i];
//    }
//}

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
	int Marker = 60; 
	clock_t start_loop, end_loop;
    for (int i = 0; i < N_particles; i++) {
        particles[i].x = gsl_rng_uniform_int(rng, Nx);
        particles[i].y = gsl_rng_uniform_int(rng, Ny);
    }
FILE *stats = fopen("stats.bin", "wb");
	if(!stats) {
		perror(" Failed to opent the binary file\n");
		free(particles); 
		return; 
	}
FILE *traj_file = fopen("trajectory.bin", "wb");
	if(!traj_file) {
		perror(" Failed to opent the binary file\n");
		free(particles); 
		return; 
	}
	for (int t = 0; t <= T_max; t++) {
	start_loop = clock(); 
        for (int i = 0; i < N_particles; i++) { 
	    int dir = gsl_rng_uniform_int(rng, 4);
            int x = particles[i].x;
            int y = particles[i].y; 	
	    double c; 
	    empirical[x][y] = empirical[x][y] + 1; 
            int xn = x, yn = y;
            if (dir == 0) {
		xn = (x + 1) % Nx;
		c = c_y[IDX(x,y)]  ; // north direction e  
		}
            else if (dir == 1) {//south
			xn = (x - 1 + Nx) % Nx;
			c = - c_y[IDX(xn, y)]; 
			}
            else if (dir == 2){//east
			 yn = (y + 1) % Ny;
			c = c_x[IDX(x,y) ]; 
			}
            else { //west
			yn = (y - 1 + Ny) % Ny;
			c = - c_x[ IDX(x, yn )];  
		}

            double V_old = V[IDX(x, y)];
            double V_new = V[IDX(xn, yn)];
            double dE = V_new - V_old;
            if (gsl_rng_uniform(rng) < (1 + c) * exp(-beta * dE)) {
                particles[i].x = xn;
                particles[i].y  = yn;
            // We save the movement with sign.  at each site[x][y], edge_x[x][y] corresponds to the edge (x,y) (x+ 1, y) and edge_y it the edge (x, y ) (x, y +1): Edi
		// EDIT: WE SAVE THE DIRECTED JUMPS FROM A SITE, THE EMPIRICAL FLOW 
		if (dir == 0)  north[x][y] =  north[x][y] + 1; //North direction 
           	else if (dir == 1)  south[x][y] = south[x][y] + 1; //south
          	else if (dir == 2) east[x][y] = east[x][y] + 1 ; //east
            	else  west[x][y] = west[x][y] +1; // west
		}
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
		if( t % snapshot_interval == 0){
			fwrite(&t, sizeof(int), 1 ,traj_file); 	
   			fwrite( particles, sizeof(Particle_State), N_particles, traj_file);
			fwrite(&Marker, sizeof(int),1, traj_file);  
			fwrite(&t, sizeof(int), 1 , stats ); 	
   			fwrite( empirical_flat, sizeof(double), Nx*Ny, stats);
   			fwrite( north_flat, sizeof(double), Nx*Ny, stats);
   			fwrite( south_flat, sizeof(double), Nx*Ny, stats);
   			fwrite( east_flat, sizeof(double), Nx*Ny, stats);
   			fwrite( west_flat, sizeof(double), Nx*Ny, stats);
			fwrite(&Marker, sizeof(int),1, stats);  
 			for( int i=0; i < Nx;  i ++){
				for(int j = 0; j < Ny ; j++){
					empirical[i][j] = 0; 
					north[i][j] = 0; 
					south[i][j] = 0; 
					east[i][j] = 0; 
					west[i][j] = 0; 
				}
			}
       			end_loop = clock();
        		double elapsed = (double)(end_loop - start_loop) / CLOCKS_PER_SEC; 
			printf(" Time taketook %.4f seconds\n", elapsed); 
   		}
	}
    fclose(traj_file);
    fclose(stats);  
}

int main() {
 //   read_config("config.txt");
    gsl_rng_env_setup();
    rng = gsl_rng_alloc(gsl_rng_mt19937);
    V = malloc(Nx * Ny * sizeof(double));
    c_x = malloc(Nx * Ny * sizeof(double));
    c_y = malloc(Nx * Ny * sizeof(double));
    mu = malloc(Nx * Ny * sizeof(double));
    double_well();
	check_orthogonality(); 
    compute_mu(); 
    gradVx = malloc(Nx * Ny * sizeof(double));
    gradVy = malloc(Nx * Ny * sizeof(double));
    empirical_flat = malloc(Nx*Ny*sizeof(double)); 
    b_east_flat = malloc(Nx*Ny*sizeof(double)); b_west_flat = malloc(Nx*Ny*sizeof(double)); 
    b_south_flat = malloc(Nx*Ny*sizeof(double)); b_north_flat = malloc(Nx*Ny*sizeof(double)); 
    north_flat = malloc(Nx*Ny*sizeof(double)); south_flat = malloc(Nx*Ny*sizeof(double)); 
    west_flat = malloc(Nx*Ny*sizeof(double)); east_flat = malloc(Nx*Ny*sizeof(double)); 
	// We allocate the infinitesimal empirical measure   
 	empirical = malloc(Nx * sizeof(double*));
	// We allocate e_1, which save the infinitesimal emprirical jumps to the right and to the north, respectively 
    east  = malloc(Nx*sizeof(double*)); 
    west = malloc(Nx*sizeof(double*)); 
    north  = malloc(Nx*sizeof(double*)); 
    south = malloc(Nx*sizeof(double*)); 
   	// Estimation of b: empirical current/ empirical measure 
	 b_east  = malloc(Nx*sizeof(double*)); 
    b_west = malloc(Nx*sizeof(double*)); 
    b_north  = malloc(Nx*sizeof(double*)); 
    b_south = malloc(Nx*sizeof(double*)); 
     for(int i = 0; i < Nx; i++ ){
	empirical[i] = &empirical_flat[ i*Ny]; 
	west[i] = &west_flat[Ny*i];
	east[i] = &east_flat[Ny*i];
	north[i] = &north_flat[Ny*i]; 
	south[i] = &south_flat[Ny*i]; 
	b_west[i] = &b_west_flat[Ny*i];
	b_east[i] = &b_east_flat[Ny*i];
	b_north[i] = &b_north_flat[Ny*i]; 
	b_south[i] = &b_south_flat[Ny*i]; 
		for (int j = 0; j < Ny; j++){
		empirical[i][j] = 0; 
		east[i][j] = 0; 
		north[i][j] = 0; 
		south[i][j] = 0;
		west[i][j] = 0;
		b_east[i][j] = 0; 
		b_north[i][j] = 0; 
		b_south[i][j] = 0;
		b_west[i][j] = 0;
		} 
	}
    mucx = malloc(Nx * Ny * sizeof(double));
    mucy = malloc(Nx * Ny * sizeof(double));
    divergence = malloc(Nx * Ny * sizeof(double));
    particles = malloc(N_particles * sizeof(Particle_State));
        
   // compute_fields();
   // write_flat_csv("V.csv", V);
   // write_flat_csv("mu.csv", mu);
   // write_flat_csv("gradVx.csv", gradVx);
   // write_flat_csv("gradVy.csv", gradVy);
   // write_flat_csv("c_x.csv", cx);
   // write_flat_csv("c_y.csv", cy);
   // write_flat_csv("muc_x.csv", mucx);
   // write_flat_csv("muc_y.csv", mucy);
   // write_flat_csv("divV_dot_c.csv", divergence);
    simulate();
    gsl_rng_free(rng);
    free(V); free(mu); free(gradVx); free(gradVy);
    free(cx); free(cy); free(mucx); free(mucy); free(divergence);
    free(particles);
    free(c_x); free(c_y); 
    free(empirical); free(empirical_flat); 
    free(east); free(west); free(north); free(south); 
    free(east_flat); free(west_flat); free(north_flat); free(south_flat); 
    free(b_east); free(b_west); free(b_north); free(b_south); 
    free(b_east_flat); free(b_west_flat); free(b_north_flat); free(b_south_flat); 
    return 0;
}

