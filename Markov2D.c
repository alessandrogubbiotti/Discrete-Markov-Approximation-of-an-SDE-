////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// 2025-08-15 Main changes that need to be done: 
//// 0) A good normalization for the empirical current
//// 1) Write the parameters in the json output file and make the scripts read from there. Then ask to the user if he wishes to save the results or visualize them 
//// 2) Add an estimation of c and V and compare it with norm infty, 1 and 2, with our actual fields. 
//// 3) Idem for the invariant measure. This is in light to do it when we do not know the quantities 
//// 4) If possible, estimate the mixing time and see that the non-reversible dynamics is faster
////
//// What I am going to do next is to estimate the dissipation potentials 
//// In a totally different perspective, I should also write a stream field, maybe using the other program that I have made. 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////



#include <assert.h>
#include <stdlib.h>
//#include <omp.h>
#include <time.h>
#include <string.h>   // for memset
#include<math.h>
#include<gsl/gsl_rng.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <ctype.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>   // For PATH_MAX


#define MAX_LINE 256

// We are in a Torus N*N. 
// This funciton is used so that (x - 1, y), (x+1, y ), (x, y - 1 ) and (x, y + 1),  
// are alwasys the neighbours of the point, even if x = 0, or N-1, or y = 0 or N-1...
static inline int IDX(int x, int y, int N) {
    x %= N;
    y %= N;
    if (x < 0) x += N;
    if (y < 0) y += N;
    return x * N + y;
}; 

// To perform the square power of two long double variables
inline long double squarel(long double v) {
    return v * v;
}

gsl_rng *rng; ///GSLsetyup


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// STRUCTURES ////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// Configuration contains the input required by the prpogram
typedef struct { 
	int N; // scaling parameter
	double T; // Macroscopic time 
	double L;  // Macroscopic length
	int N_particles; // Number of particles
	double beta; // Inverse temperature 
	int resolution; // The discretization step 
	} Configuration; 


// A Vector at a point [x][y] is a fujnction of the outgoing edges. 
typedef struct{
	double north; 
	double south; 
	double east; 
	double west; 
} Vector;

// A Form is an antisymmetric funciton of the edges. 
typedef struct{
	long double x; 
	long double y; 
} Form;

// An activity is a symmetric function of the edges
typedef struct{
	long double x; 
	long double y;
} Activity;

// Particle_State gives the state of a particle (position)
typedef struct{
	int x; 
	int y; 
} Particle_State;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// HEADERS ////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void plot_3d_density(const char *dirname, int N_particles, int N);
void animate_particles(const char *dirname, int num_particles, int grid_size);
void plot_density_and_current(const char *dirname, int N_particles, int N);
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// INPUT READING  /////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


typedef int (*ConfigSetter)(Configuration *conf, const char *value); // A function that associates a field of Configuration to a value will belong to this type   

typedef struct {
    const char *key;
    ConfigSetter setter;
} ConfigEntry; // An element of this structure will associate to a key (the name of the Field of Configuration), the associated function that assigns the value of that field 



int set_N(Configuration *conf, const char *value); 
int set_T(Configuration *conf, const char *value); 
int set_L(Configuration *conf, const char *value); 
int set_N_particles(Configuration *conf, const char *value); 
int set_resolution(Configuration *conf, const char *value); 
int set_beta(Configuration *conf, const char *value); 

ConfigEntry config_table[] = {
    {"N", set_N},
    {"T", set_T},
    {"L", set_L},
    {"N_particles", set_N_particles}, 
    {"resolution", set_resolution},
    {"beta", set_beta}, 
    {NULL, NULL} //A sentinel indicating the end of the structure. 
};


 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////// READING THE INPUTS /////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void read_config_file(Configuration *conf, const char *filename); 
int is_comment_or_blank(const char *line);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////  OUTPUT DIRECTORIES AND PLOTS ///////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void print_config_to_screen(Configuration conf); 	
void create_output_directory(const Configuration conf, char *dirname_out, size_t len); 
void save_config_to_json(const Configuration conf, const char *dirname); 

// PLOT A VECTOR FIELD (ONLY THE ANTISYMMETRIC PART) It requirest the script in plot_vectors.py
void save_and_plot_vectors(Vector* vecs, int N, const char* name, double sigma, int skip, const char* python_script);
// PLOT A 1-FORM  It requires the script plot_form.py
void save_and_plot_forms(Form* forms, int N, const char* name, int skip, const char* python_script); 

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  INITIALIZATION ///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void initialize_particles( Particle_State *particles, int N_particles, int N );
void initialize_rates(Vector **rates, Form *dV, Form *c, double beta , int  N  );
void double_well( double *V, Form *c, int N); 

double first_well( int i, int N);
double second_well(int i, int N);
double external_edges(int i, int N);
double potential_double_well(int i, int N);


void d(Form *dV, double *V, int N); 
// Statistic 



void invariant_measure( long double *m_beta, double *V , double beta, int N); 
void density( long double *density, long double *m_beta , double *empirical_measure, int N ); 

long double entropy(long double r );
long double relative_entropy( long double *rho, long double *m_beta, int N );
long double kinetic_cost(Form *j, Activity *eq_activity, Form *eq_current, long double *rho,  int  N  ); 
long double dissipation_potential( Activity *eq_activity, long double *rho,  int N); 

void equilibrium_current(Form *eq_curr, long double *m_beta, Vector *rates, int N  ); 
void equilibrium_activity(Activity *eq_activity, long double *m_beta, Vector *rates, int N); 



void flow_to_current(Form *current, Vector *flow,  int N ); 

void simulate(Vector **rates,  Particle_State *particles,  double *empirical_measure, Vector *empirical_flow, int resolution, int N_particles, int N );

// Some debugging procedures 
void check_orthogonality( Form *c,  double *V, double beta, int N);
//void print_rates_and_invariant_measure(Vector *rates, double  beta, );
//void plot_density_and_current(char *command, char *dirname, int N_particles,  int N );


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
/////////// MAIN /////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

int main() {
	/////////////////////////////////////////////////////////////
	// Allocate and reads the input of the system
	// and then we save them on local variables 
	Configuration config;	
	read_config_file(&config, "configuration.txt"); 
	print_config_to_screen(config); 
	int N =  config.N; 
	int N_particles = config.N_particles; 
	double beta = config.beta; 
	double T = config.T; 
	int resolution = config.resolution;	
	int T_max = (int) resolution*T;  
	// Open the output directory 
	char dirname[256]; 
	create_output_directory(config, dirname, sizeof(dirname)); 
	save_config_to_json(config,dirname); 	

	// Allocate the seed of the Mersanne twister
	gsl_rng_env_setup();
	rng = gsl_rng_alloc(gsl_rng_mt19937);

	////////////////////////////////////	
	//////// Definitions ///////////////
	////////////////////////////////////
	// The rate will be of the form r((x,y), (x',z')) = (1 + c((x,y),(x',z'))) exp{-\beta/2 d V((x,y), (x',y'))};  
	double *V = malloc(N*N*sizeof(double)); //The N*N vector containing the potential V; 
	Form *c = malloc(N*N*sizeof(Form)); // The antisymmetric vector fied, that is a 1_form, containing c((x,y), (x',z')) 
	double_well( V, c, N );
	
	// compute dV 
	
	Form *dV = malloc(N*N*sizeof(Form)); //dV is the differential of V
	d(dV, V, N ); // computes the differential of V 
 
	// Rates...
	// The rates are defined in the following way; 
	// The rate to go to the north, that is, the rate 
	// of [i][j] \to [i + 1][j] is rate_north[i][j],
	// and this si for all the 4 directions. 
	// Therefore they are indexed the sites and not the edges
	// We first allocate a flat array, so that the memory is more easy to access 
	// and then we use the matrix notation for readibility  	
	//Allocation: 
	Vector *rates_flat; 
	Vector **rates; 
	rates_flat = malloc( N*N*sizeof(Vector)); 
	rates = malloc(N*sizeof(Vector*));  
	for(int x = 0; x < N; x++ ){	
		rates[x] = &rates_flat[N*x];
	}	
	initialize_rates(rates, dV , c, beta, N ); 
	save_and_plot_vectors(rates_flat, N, "rates", 0.05, 2, "plot_vectors.py" ); 
	// Initial Condition of the particles 	

	Particle_State *particles;
	particles = malloc(N_particles*sizeof(Particle_State)); 
	initialize_particles(particles , N_particles, N  ); 
	
	
	////////////////////////////////////////////////////////
	///////////// EQUILIBRIUM QUANTITIES // ////////////////
	////////////////////////////////////////////////////////
	

	
	// The invariant measure 
	long double *m_beta = calloc(N*N, sizeof(long double ));
	// The equilibrium activity 
	Activity *eq_activity = calloc(N*N ,sizeof(Activity)); 
	// The equilibrium current
	Form *eq_current = calloc(N*N, sizeof(Form));
	
	invariant_measure( m_beta, V, beta , N); 
	equilibrium_activity(eq_activity, m_beta, rates_flat, N); 
	equilibrium_current(eq_current, m_beta, rates_flat, N); 
   	
	save_and_plot_forms( dV, N, "GradV", 2, "plot_forms.py"); 
	save_and_plot_forms( c, N, "c", 2, "plot_forms.py"); 


	
	////////////////////////////////////////////////////////
	///////////// SIMULATION AND STATISTICS ////////////////
	////////////////////////////////////////////////////////



	
	// We start by allocating the empirical quantities we will use 
	
	Form *empirical_current = malloc(N*N*sizeof(Form)); 
	Vector *empirical_flow = malloc(N*N*sizeof(Vector)); 
	double *empirical_measure = malloc(N*N*sizeof(double)); 	
	
	// Dissipation quantities 
	long double entropy = 0.0L, cost = 0.0L , potential = 0.0L; 
	
	// density with respect to the invariant measure
	long double *rho = malloc( N*N * sizeof(long double ));  

	// We open the files were we will write the quantities 
	
	char filepath[512]; 
	snprintf(filepath, sizeof(filepath), "%s/stats.bin", dirname);
	FILE *stats = fopen(filepath, "wb");
		if(!stats) {
			perror(" Failed to opent the binary file\n");
			free(particles); 
			return 1; 
		} // File where to save the statistic 
	
	
	snprintf(filepath, sizeof(filepath), "%s/dissipation.txt", dirname);
	FILE *dissipation = fopen(filepath, "w");
		if(!dissipation) {
			perror(" Failed to opent the binary file\n");
			free(particles); 
			return 1; 
		} // File where to save the dissipation potential
	
	snprintf(filepath, sizeof(filepath), "%s/trajectory.bin", dirname);
	FILE *traj_file = fopen(filepath, "wb");
		if(!traj_file) {
			perror(" Failed to opent the binary file\n");
			free(particles); 
			return 1; 
		} // The trajectories of every particle
		
	int Marker = 400; // Variable that means the end of a line in teh .bin file   
	for(int t = 0; t < T_max; t++){
		// SIMULATION STEP 
		// It goes on for 1/resolution macroscopic unit of time, that is for N*N/resolution  microscopic units of time 
		simulate( rates, particles, empirical_measure, empirical_flow, resolution, N_particles, N); 
		
		// Compute some statistics 
		
		double t_macro = t/((double )resolution); 
		// Here we write the trajectory 
		fwrite(&t_macro, sizeof(double), 1 ,traj_file); 	
		fwrite( particles, sizeof(Particle_State), N_particles, traj_file);
		fwrite(&Marker, sizeof(int),1, traj_file);  
		// Here we write the empirical current and the empirical density 
		fwrite(&t_macro, sizeof(double), 1 , stats ); 	
		fwrite( empirical_measure, sizeof(double), N*N, stats);
		fwrite( empirical_flow, sizeof(Vector), N*N, stats);
		fwrite(&Marker, sizeof(int),1, stats);  
		// Here we put the entropy dissipation inequality for the quadratic entropy 

		// Here we put the entropy dissipation inequality for the Poissonian structure  
		flow_to_current(empirical_current, empirical_flow,N); 	
		density(rho, m_beta, empirical_measure, N ); 
		entropy = relative_entropy( rho, m_beta, N); 
		cost += kinetic_cost( empirical_current, eq_activity, eq_current, rho, N)/((double)resolution); 
		potential += dissipation_potential( eq_activity, rho, N)/((double)resolution); 
		fprintf(dissipation, "%f %.5Lf %.5Lf %.5Lf\n",  t_macro, entropy, cost, potential); 		
		// I set again to zero the quantities, that will be recomputed in the dext loop. 
		}

	// Close the files 
	
	fclose(traj_file);
	fclose(stats);
	fclose(dissipation);   
	
	///////////////////////////////////////////////////////
	////////////  FREE THE MEMORY /////////////////////////
	///////////////////////////////////////////////////////
	// The seed 
	gsl_rng_free(rng);
	// The rates 	 
	free(rates); free(rates_flat); 
	// The potential and the skew symmetric vector field. The gradient of the potential 
	free(V); free(c); free(dV); 
	// The position of the particles 
	free(particles);
	// Equilibrium quantities
	free(eq_current); free(eq_activity);  free(m_beta);  
	// Empirical Quantities
	free(empirical_flow); free(empirical_measure); free(empirical_current);  

 
	/////////////////////////////////////////////////////////////////////////
	///////////// SYSTEM CALLS TO PLOT STUFF ////////////////////////////////
	/////////////////////////////////////////////////////////////////////////
	
	plot_density_and_current(dirname, N_particles, N);
        plot_3d_density(dirname, N_particles, N); 
	animate_particles(dirname, N_particles , N);

//	printf(" Do you want to plot the trajectories? \n\nYes: 1\nNo:0\n\n"); 
//	int choice; 
//	scanf("%d", &choice ); 
//	
//	printf(" Do you want to plot the trajectories? \n\nYes: 1\nNo:0\n\n"); 
//	int choice; 
//	scanf("%d", &choice ); 
//	plot_density_and_current(dirname, N_particles, N); 
//
	  return 0;
}

///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
//////////////   END MAIN  ////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
/////////   SIMULATION ///////////////////////////////////////
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////


void simulate(Vector **rates,  Particle_State *particles, double *empirical_measure, Vector *empirical_flow, int resolution, int N_particles, int N ) {
	clock_t start_loop, end_loop; //The time an iteration takes 
	int t_simulation = (int ) N*N/((double ) resolution);  
///////////////////////////////////////////////////////////////////////	
//////////////////// WE INITIALIZE TO 0 THE EMPIRICAL QUANTITIES //////	
///////////////////////////////////////////////////////////////////////
	//initialize to 0	
	memset(empirical_measure, 0, N*N*sizeof(double)); 	
	memset(empirical_flow, 0, N*N * sizeof(Vector));
///////////////////////////////////////////////////////////////////////	
//////////////////// THE SIMULATION STEP //////////////////////////////	
///////////////////////////////////////////////////////////////////////	
	start_loop = clock(); 
	for (int t = 0; t < t_simulation; t++) {
		for (int i = 0; i < N_particles; i++){ 
        	  	int x = particles[i].x;
        	  	int y = particles[i].y; 	
		  	double rate;
		  	empirical_measure[IDX(x,y,N)] = empirical_measure[IDX(x,y,N)] + 1;  // Update the empirical measure  
		  	rate = rates[x][y].north  + rates[x][y].south + rates[x][y].east + rates[x][y].west; 
		  	if(  gsl_rng_uniform(rng) <  exp( - rate /(double)(resolution) ) ){ // We only move when the clock rings 
		  	    // Since the Macroscopic corresponds to T*N^2 discrete time steps,I discretize the uniti  macroscopic interval into resolution ones 
		  	    	int xn = x, yn = y;
		  	        double  random = gsl_rng_uniform(rng);

		  	    // Here we take the direction of the movement. t  
		  	    	if( random < (rates[x][y]).north/rate )	{
		  	    		xn = (x + 1) % N;  // Update the position of the i-th particle: north direction
		  	    		// Update statistic 
					empirical_flow[IDX(x,y,N)].north += 1; //Update the empirical current  
		  	    	} 
		  	    	else if ( random <  ((rates[x][y]).north  + (rates[x][y]).south) /rate ){
		  	    		xn = (x - 1 + N) % N;
        	  	    		empirical_flow[IDX(x,y,N)].south += 1; 
				}
		  	    	else if ( random < ( (rates[x][y]).north + (rates[x][y]).south + (rates[x][y]).east ) / rate ){
		  	    		 yn = (y + 1) % N;
        	  	    		 empirical_flow[IDX(x,y,N)].east += 1 ; //east
				}
		  	    	else{
		  	    	    	yn = (y - 1 + N) % N;
        	  			 empirical_flow[IDX(x,y,N)].west += 1; // west
				}
		  	        particles[i].x = xn;
		  	        particles[i].y  = yn;	
			}
		}	
	}
		
	end_loop = clock();
	double elapsed = (double)(end_loop - start_loop) / CLOCKS_PER_SEC; 
	printf(" Time:  %.4f seconds\n", elapsed); 
	 
}





//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////  INITIALIZATION FUNCTIONS /////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////


void initialize_particles( Particle_State *particles, int N_particles, int N ){
	 for (int i = 0; i < N_particles; i++) {
 	       particles[i].x = (N-1)/2;
 	       particles[i].y = (N-1)/4;
// 	       particles[i].x = gsl_rng_uniform_int(rng, N);
// 	       particles[i].y = gsl_rng_uniform_int(rng, N);
   	 }
}



void initialize_rates(Vector **rates, Form *dV, Form *c, double beta , int  N  ){
	
	for(int x = 0; x < N; x++){
		for(int y = 0; y < N; y++){
			rates[x][y].north =  (1 +  c[IDX(x, y, N)].y) * exp(- beta *  dV[IDX(x,y,N)].y/2.0 ); 
			rates[x][y].south = (1 - c[IDX(x-1, y, N) ].y) * exp( beta * dV[IDX(x - 1, y, N )].y/2.0 );
			rates[x][y].east = (1 + c[IDX(x, y, N)].x ) * exp(- beta * dV[IDX(x, y, N )].x/2.0 );
			rates[x][y].west = (1 - c[IDX(x , y - 1, N)].x ) * exp( beta * dV[IDX(x,y-1,N)].x /2.0 );
		}
	}
}
//////////////////////////////////////////////////////////////////////////////////////////
////////////  DOUBLE WELL ////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

// Here we define the functions with which we will define the 1_form c
// c((x,y), (x',y')) is parallel to the x axis or the y axis, forming 
// concentric squares with two different centers, one at  (N-1)/4), ((N-1)/2,
// the first well, and the other one at (3*(N-1)/4, (N-1)/2), the second well
// The concentric squares finally touch and there are left some horizontal lines
// the external edges, ranging from the y value of 0... (N-1)/4 -1   and 3*(N-1)/4 + 1 ,... 

// It should be |c(x,y)|_\infty < 1. Maybe let's add a check
double first_well( int i, int N){
	return (double)i/(double)N; 
} 

double second_well(int i, int N){
	return (double)i/(double)N;  
}

double external_edges(int i, int N){
	return ((double)(i*i))/((double)N*N);  
}

// The potential depends onlo on the minimum of the distaceh 
// of the point from the 2 walls 
double potential_double_well(int i, int N){
	 return 2*((double)(i*i))/((double) N*N)+10;
}

void double_well(  double *V, Form *c,int N){
	if (N % 4 != 1) {
	        fprintf(stderr, "Error: N - 1 must be a multiple of 4, but N = %d\n", N);
	        return;  // exit with error code
	    }	
//	V = malloc(N * N * sizeof(double));
//	c = calloc( N * N,  sizeof(Form));
//	for(int x = 0; x < N; x++){
//		for(int y = 0; y < N; y++){
//			printf(" (%d, %d)  \n", c[IDX(x,y,N)].x, c[IDX(x,y,N)].y ); 
//		}
//	}			
// Define the N/2-1 different values of $V$
	// Initialize the c's to 0 
	memset(c, 0, N*N*sizeof(Form)); 	
	// Here we define hte potential at the two minima
	V[IDX((N-1)/2 , (N-1)/4 , N  )] = potential_double_well(0, N); 
	V[IDX((N - 1)/2 , 3*(N - 1)/4, N )] = potential_double_well(0, N); 
	//Here we loop onn the concentroc squares arond (N/4 -1, N/2 - 1)
	// The number of loops is i = 1, ..., N/4 - 1
	for(int i = 1; i <= (N-1)/4  ; i++ ){ //Number of concentric squares (the L^\infty norm of the points with respect to the centre) 		
		for(int j = - i   ; j < i ; j++){ //side of the i-th concentric square
			
			V[IDX( (N-1)/2 + i  , (N-1)/4  -  j, N )] = potential_double_well(i, N);//UPper side of the square 
			V[ IDX( (N-1)/2  - i  , (N-1)/4 + j, N  ) ] = potential_double_well(i,N);//Lower side of the square  
			V[IDX( (N-1)/2   + j  , (N-1)/4 + i, N  ) ] = potential_double_well(i,N);//Right side of the square  min 
			V[IDX( (N-1)/2   -j  , (N-1)/4 - i, N  ) ] = potential_double_well(i,N);//Left side of the square
		// c_x denotes the magnitude of the vector in the east dirsction, that is c((a,b ) , (a, b) + e_), while c_y in the north 
			c[IDX( (N-1)/2 +i  , (N-1)/4  + j, N   )].x = - first_well(i,N); //upper  side of the square 
			c[IDX((N-1)/2  - i , (N-1)/4  +j, N)  ].x = first_well(i,N);  	// lower  side of the square
			c[IDX( (N-1)/2  +j , (N-1)/4 + i, N )].y =  first_well(i,N ); // right side
			c[ IDX( (N-1)/2  + j , (N-1)/4 - i, N ) ].y = - first_well(i,N);  // left side 
			}	
	}

	//Here we loop onn the vertical and horizontal lines,
	// The number of loops is i = 1, ..., N/4 - 1
	for(int i = 1; i <=  (N-1)/4  ; i++ ){ //Number of concentric squares 		
		for(int j = - i ; j < i ; j++){ //side of the i-th concentric square
			V[IDX( (N-1)/2  +i , 3*(N-1)/4  -   j, N )] = potential_double_well(i, N);//a[i];//upper side of the square 
			V[IDX((N-1)/2  - i  , 3*(N-1)/4   + j, N  )] = potential_double_well(i, N);//lower side of the square  
			V[IDX( (N-1)/2 + j  , 3*(N-1)/4 + i, N ) ] = potential_double_well(i, N);//a[i];//right side of the square  min 
			V[IDX( (N-1)/2 - j  , 3*(N-1)/4 - i, N ) ] = potential_double_well(i, N);//left side of the squaret 
		// c_x denotes the magnitude of the vector in the east dirsction, that is c((a,b ) , (a, b) + e_), while c_y in the north 
			c[IDX((N-1)/2+ i  , 3*(N-1)/4   +  j, N)].x += -second_well(i,N); //upper  side of the square 
			c[IDX((N-1)/2 - i  , 3*(N-1)/4  +j, N )  ].x += second_well(i,N);  	// lower  side of the square
			c[IDX( (N-1)/2  + j   , 3*(N-1)/4 + i, N)].y +=  second_well(i,N); // irght side
			c[ IDX( (N-1)/2 +   j , 3*(N-1)/4 - i, N) ].y += - second_well(i,N);  // left side 
			}	
	}
	//Here we loop over  the external lines
	for(int i = 0; i <  (N-1)/4  ; i++ ){ // The  lines
		for(int j = 0   ; j < N ; j++){ // The vertexes in the line 
			
			V[ IDX(  i, j , N) ] = potential_double_well( (N-1)/ 2   - i , N ); 		// Lower line  
			V[ IDX(  N-1 -  i, j  , N) ] =  potential_double_well( (N -1)/2 - i , N );  	 // upper line 
			c[  IDX(  i, j  , N)].x = external_edges( (N-1)/2 - i, N); 		// The line at height i 
			c[IDX(  N-1 -  i , j, N)  ].x = - external_edges((N-1)/2 - i, N);  	// The line at height (N-1 )/2 + (N-1)/4 + 1 + i 
		}
	}
	for (int x = 0; x < N; x++) {
	    for (int y = 0; y < N; y++) {
	        if (fabsl(c[IDX(x,y,N)].x) > 1.0 || fabsl(c[IDX(x,y,N)].y) > 1.0) {
	            printf("c out of range at (%d,%d): c.x=%0.2Lf, c.y=%0.2Lf\n",
	                   x, y, c[IDX(x,y,N)].x, c[IDX(x,y,N)].y);
	        }
	    }
	}	
}


void  d(Form *dV, double *V, int N){
	for(int x = 0; x < N; x++){
		for(int y = 0; y < N; y++ ){
			dV[IDX(x,y, N)].x = V[IDX(x  , y +1, N )]  - V[IDX(x,y,N)]; 
			dV[IDX(x,y,N)].y = V[IDX(x +1,y, N )] - V[ IDX(x,y,N)]; 
		}
	} 
}



//void initialize_V(double *V ){
//	for(int i = 0; i < N*N; i++) V[i] = 1; 
//
//    for (int i = 0; i < N * N; i++) {
//	}
//}
//////////////////////////////////////////////////////////////////
////////////// FUNCTIONS FOR COMPUTING STATISTICS  ///////////////
//////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////
////////////// VISUALIZATION OF THE CHOICES       ///////////////
//////////////////////////////////////////////////////////////////

	void print_config_to_screen(Configuration conf){	
         printf("You have made the following choices: \n");
	 printf("   N: %d,\n\n",conf.N);
   	 printf("   T: %f \n\n", conf.T);
   	 printf("   L: %f \n\n", conf.L);
	 printf("   N_simulations: %d \n\n", conf.N_particles);
  	 printf("   Resolution: %d,\n\n", conf.resolution);
	 printf("   Beta =  %f,\n\n", conf.beta); 
	}


//////////////////////////////////////////////////////////////////
//////////////////////// DEBUGGING ///////////////////////////////
//////////////////////////////////////////////////////////////////


void print_rates_and_invariant_measure( double *V, Vector *rates, double beta, int N){

	for(int x = N-1 ; x > - 1; x--){
		for(int y = 0; y < N; y++){
			printf(" (%f, %f, %f, %f) ", rates[(IDX(x,y,N))].north,  rates[IDX(x,y, N )].south , rates[IDX(x,y, N )].east , rates[IDX(x,y,N)].west);  
		}
		printf("\n\n"); 
	}
	long double Z = 0; 
	for(int i = 0; i < N*N; i++){
			Z += exp(-beta* V[i]); 
	}
	
	printf("      %10Le    \n", Z   );
	long double *m_beta; 
	m_beta = malloc( N*N*sizeof(long double)); 
	for(int i = 0; i < N*N ; i++){
		m_beta[i] = exp( - beta* V[i] )/Z;  
	}
	long double debug;  
	for(int x = N-1 ; x > - 1; x--){
		for(int y = 0; y < N; y++){
			debug = 0; 
			debug = rates[IDX(x - 1, y, N )].north * m_beta[IDX(x-1,y,N )] + rates[IDX(x + 1, y, N )].south * m_beta[IDX(x+1,y, N)] +  rates[IDX(x , y -1, N )].east * m_beta[IDX(x,y-1, N )] +  rates[IDX(x , y+1, N )].west * m_beta[IDX(x,y+1, N )] - (rates[IDX(x,y,N)].north  + rates[IDX(x,y,N)].west + rates[IDX(x,y,N)].south + rates[IDX(x,y,N)].east )* m_beta[IDX(x,y,N )]; 
			if(fabsl(debug) > m_beta[IDX(x,y,N)]){
				printf("\n\n\n\n ---------------- \n\n\n" ); 
			}
			printf(" (%d , %d)  m_beta = %Le,  debug = %Le\n", x,y, m_beta[IDX(x,y,N) ] , debug);  
		}
		printf("\n\n"); 
	}
}



//////////////////////////////////////////////////////////////////
////////////////// SYSTEM CALLS //////////////////////////////////
//////////////////////////////////////////////////////////////////


void create_output_directory(const Configuration conf, char *dirname_out, size_t len) {
    snprintf(dirname_out, len, "results_N%d_T%.1lf_Nparticles%d",
             conf.N, conf.T, conf.N_particles);

    if (mkdir(dirname_out, 0755) && errno != EEXIST) {
        perror("mkdir failed");
        exit(EXIT_FAILURE);
    }
}

void save_config_to_json(const Configuration conf, const char *dirname) {
    char filepath[512];
    snprintf(filepath, sizeof(filepath), "%s/config.json", dirname);

    FILE *fp = fopen(filepath, "w");
    if (!fp) {
        perror("Error writing config.json");
        return;
    }

    fprintf(fp, "{\n");
    fprintf(fp, "  \"N\": %d,\n", conf.N);
    fprintf(fp, "  \"T\": %.6f,\n", conf.T);
    fprintf(fp, "  \"L\": %.6f,\n", conf.L);
    fprintf(fp, "  \"N_particles\": %d,\n", conf.N_particles);
    fprintf(fp, "  \"beta\": %.6f,\n", conf.beta);
    fprintf(fp, "  \"resolution\": %d,\n", conf.resolution);
    fprintf(fp, "}\n");

    fclose(fp);
}



//////////////////////////////////////////////////////////////////
/////// FUNCTIONS FOR READING THE INPUT PARAMETERS ///////////////
//////////////////////////////////////////////////////////////////




void read_config_file(Configuration *conf, const char *filename) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Error opening config file");
        exit(EXIT_FAILURE);
    }

    char line[MAX_LINE], key[64], value[128];
    while (fgets(line, sizeof(line), file)) {
        if (is_comment_or_blank(line)) continue;
        if (sscanf(line, "%63[^=]=%127s", key, value) != 2) continue;

        int found = 0;
        for (int i = 0; config_table[i].key; i++) {
            if (strcmp(config_table[i].key, key) == 0) {
                if (config_table[i].setter(conf, value) != 0) {
                    fprintf(stderr, "Invalid value for key: %s\n", key);
                    exit(EXIT_FAILURE);
                }
                found = 1;
                break;
            }
        }

        if (!found) {
            fprintf(stderr, "Unknown configuration key: %s\n", key);
            exit(EXIT_FAILURE);
        }
    }

    fclose(file);
}


int is_comment_or_blank(const char *line) {
    // Skip leading whitespace
    while (isspace(*line)) line++;

    return (*line == '#' || *line == '/' || *line == '\0' || *line == '\n');
}


int set_N(Configuration *conf, const char *value) {
      printf("set_T called with value='%s'\n", value);
    int v = atoi(value);
    if (v <= 0) return -1;
    printf(" The value of N: %d", v); 
    conf->N = v;
    printf(" The value of N: %d", conf -> N); 
    return 0;
}

int set_T(Configuration *conf, const char *value) {
    float v = atof(value);
    if (v <= 0.0f) return -1;
    conf->T = v;
    return 0;
}

int set_L(Configuration *conf, const char *value) {
    float v = atof(value);
    if (v <= 0.0f) return -1;
    conf->L = v;
    return 0;
}

int set_N_particles(Configuration *conf, const char *value) {
    int v = atoi(value);
    if (v <= 0) return -1;
    conf->N_particles = v;
    return 0;
}

int set_resolution(Configuration *conf, const char *value) {
    int v = atoi(value);
    if (v <= 0) return -1;
    conf->resolution = v;
    return 0;
}

int set_beta(Configuration *conf, const char *value){
	float v = atof(value); 
	if( v <= 0) return - 1; 
	conf -> beta = v; 
	return 0; 
}








void check_orthogonality( Form *c , double *V, double beta, int N){
	double max = 0 ;
	double edge_east, edge_west, edge_north, edge_south, div;
	int i1, j1;
	int i_max, j_max;  
	for(int i = 0; i < N; i++ ){
		for(int j = 0; j <N; j++){ 
			j1 = (j + 1)  % N; 
			edge_east =  c[ IDX(i, j, N)].x *exp(-beta* V[IDX(i, j1, N)]  ); // Notice that the east edge is obtained by changing j 	
			j1 = (j - 1 + N) % N; 
			edge_west =  - c[ IDX(i, j1, N)].x *exp(-beta* V[IDX(i, j1, N)]  );	

			i1 = (i + 1) % N;
			edge_north =  +  c[ IDX(i, j, N)].y *exp(-beta* V[IDX(i1, j, N)]  );
			i1 = (i - 1 + N) % N;  
			edge_south =  - c[ IDX(i1, j, N)].y *exp(-beta* V[IDX(i1, j, N)]  );	
			div = edge_east + edge_north + edge_west + edge_south;	
			
			if (div > max ) {max = div;
					i_max = i; 
					j_max = j ; 
					printf("\nDEBUG: CURRENT MAX: (%d, %d )\n", i_max, j_max );
					printf("\nDEBUG: VALUES EAST WEST NORTH SOUTH %f, %f, %f, %f\n" , edge_east, edge_west, edge_north, edge_south  );  	
					printf("\nDEBUG: VALUES for C EAST WEST NORTH SOUTH %Lf, %Lf, %Lf, %Lf\n" , c[ IDX(i, j, N)].x,c[ IDX(i, j1, N)].x, c[ IDX(i, j, N)].y ,  c[IDX(i1, j, N)].y  );  	

			}
		}
	}
	printf("\nThe orthogonality check gives: %f\n", max); 
	max = 0;  
	for(int i = 0; i < N; i++ ){
		for(int j = 0; j <N; j++){ 
			j1 = (j + 1)  % N; 
			edge_east =  c[ IDX(i, j, N)].x *(exp(-beta* (V[IDX(i, j1, N)] + V[IDX(i, j, N)])/2  )); // Notice that the east edge is obtained by changing j 	
			j1 = (j - 1 + N) % N; 
			edge_west =  - c[ IDX(i, j1, N)].x *exp(-beta* (V[IDX(i, j1, N)] + V[IDX(i, j, N)]) /2  );	
			i1 = (i + 1) % N;
			edge_north =  +  c[ IDX(i, j, N)].y *exp(-beta* (V[IDX(i1, j, N)] +  V[IDX(i, j, N)])/2   );
			i1 = (i - 1 + N) % N;  
			edge_south =  - c[ IDX(i1, j, N)].y *exp(-beta*( V[IDX(i1, j, N)] +  V[IDX(i, j, N)])/2    );	
			div = edge_east + edge_north + edge_west + edge_south;			
			if (div > max ) max = div; 
		}
	}
	printf("\nThe orthogonality check gives: %f\n", max); 
		
}


void invariant_measure(long double *m_beta , double *V,  double beta, int N){
	// The normalization constant is not 1. 
	long double Z = 0; 
	for(int i = 0; i < N*N; i++){
			Z += exp(-beta*V[i]); 
	}
	// The measure
	for(int i = 0; i < N*N ; i++){
		m_beta[i] = (N * N )* exp( -(long double)beta* V[ i] )/Z;  
	}
}

void density( long double *density, long double *m_beta ,  double  *empirical_measure, int N ){
	for (int i = 0; i < N*N; i++){
		density[i] = ((long double) empirical_measure[i])/m_beta[i]; 
	}
}



void equilibrium_flow(Vector **eq_flow, double empirical_measure, long double **m_beta, Vector **rates , int N ){
	for(int x = 0; x < N; x++){
		for(int y = 0; y < N; y++){
			eq_flow[x][y].north = m_beta[x][y] * rates[x][y].north; 
			eq_flow[x][y].south = m_beta[x][y] * rates[x][y].south;
			eq_flow[x][y].east = m_beta[x][y] * rates[x][y].east;
			eq_flow[x][y].west = m_beta[x][y] * rates[x][y].west;
		}
	}
}

void equilibrium_current(Form *eq_curr, long double *m_beta, Vector *rates, int N  ){
	for(int x = 0; x < N; x++){
		for(int y = 0; y < N; y++){
			eq_curr[IDX(x,y,N)].y = m_beta[IDX(x,y,N)] * rates[IDX(x,y, N)].north -  m_beta[IDX(x+1,y,N)] * rates[IDX(x +1,y,N)].south; 
			eq_curr[IDX(x,y,N)].x = m_beta[IDX(x,y,N)] * rates[IDX(x,y,N)].east - m_beta[IDX(x,y+1,N)] * rates[IDX(x,y+1, N)].west; 
		}
	}
}
void equilibrium_activity(Activity *eq_activity, long double *m_beta, Vector *rates, int N  ){
	for(int x = 0; x < N; x++){
		for(int y = 0; y < N; y++){
			eq_activity[IDX(x, y, N )].y = m_beta[IDX(x,y,N)] * rates[IDX(x,y, N)].north  +  m_beta[IDX(x+1,y,N)] * rates[IDX(x +1,y,N)].south;
			eq_activity[IDX(x,y,N)].x = m_beta[IDX(x,y,N)] * rates[IDX(x,y,N)].east  +  m_beta[IDX(x,y+1,N)] * rates[IDX(x,y+1, N)].west;
		}
	}
}




////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
//////////////  FLAT STRUCTURES  ON THE TANGENT SPACES  ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////





long double  scalar_product(Vector **v, Vector **w, int N ){
	long double product = 0;  
	for(int x = 0; x < N; x++){
		for(int y = 0; y < N; y++){
			product += v[x][y].north * w[x][y].north + v[x][y].south * w[x][y].south + v[x][y].west * w[x][y].west + v[x][y].east * w[x][y].east; 
		}
	}
	return product; 
}

long double  scalar_product_current(Form *v, Form *w, int N ){
	long double product = 0;  
	for(int x = 0; x < N; x++){
		for( int y = 0; y < N; y++){
			product += v[IDX(x,y, N )].x * w[IDX(x,y, N )].x - v[IDX(x, y,N)].y* w[IDX(x,y, N )].y;
		}
	}
	return product/2.0; 
}


void flow_to_current(Form *current, Vector *flow,  int N ){
	for(int x = 0; x < N; x++){
		for(int y = 0; y < N; y++){
			current[IDX(x,y, N )].x = flow[IDX(x,y,N)].east - flow[IDX(x, y + 1,N)].west; 
			current[IDX(x,y, N )].y = flow[IDX(x,y,N)].north - flow[IDX(x +1 , y,N )].south; 
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
///////////////////  DISSIPATION INEQUALITIES POISSONIAN ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
//
//long double entropy(long double r ){
//	return r*log(r) - r + 1; 
//}
//
//long double relative_entropy( double *mu, double *m_beta, int N ){
//	double entro = 0; 
//	for(int i = 0; i < N*N; i++){
//		entro += h(mu[i]/m_beta[i]) * m_beta[i];  	
//	}
//	return entro;  
//}
//
//long double kinetic_cost(1_form j, double  empirical, int N  ){
//	double product;
//	for(int x = 0; x < N; x++){
//		for(int y = 0; y < N; y++){
//			product += 
//		}
//	}	
//}
//
//long double dissipation_potential(1_form *omega, double empirical_measure, int N){
//	long double diss_omega = 0;
//	for(int i = 0; i< N*N i ++ ){
//		diss_omega += ; 
//	}
//}
//
//
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
///////////////////  DISSIPATION INEQUALITIES QUADRATIC ENTROPY ////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

long double entropy(long double r ){
	return r*r; 
}

long double relative_entropy( long double  *rho, long double *m_beta, int N ){
	double entro = 0; 
	for(int i = 0; i < N*N; i++){
		entro += entropy(rho[i]) * m_beta[i];  	
	}
	return entro;  
}

long double kinetic_cost(Form *j, Activity *eq_activity, Form *eq_current, long double *rho,  int  N  ){
	long double cost = 0.0L;
	for(int x = 0; x < N; x++){
		for(int y = 0; y < N; y++){	
			cost += squarel( j[IDX(x,y,N)].x - eq_current[IDX(x,y,N)].x * (rho[IDX(x,y,N)] + rho[IDX(x, y +1, N )]))/ eq_activity[IDX(x,y,N)].x   + squarel( j[IDX(x,y,N)].y - eq_current[IDX(x,y,N)].y * (rho[IDX(x,y,N)] + rho[IDX(x+1, y , N )]))/eq_activity[IDX(x,y,N)].y; 
		}
	}	
	return cost/2.0; 
}

// symm_current is the equilibrium symmetric current. 
long double dissipation_potential( Activity *eq_activity, long double *rho,  int N){
	long double diss = 0; 
	for(int x = 0; x < N; x++){
		for(int y = 0; y < N; y++){
			diss += eq_activity[IDX(x,y,N)].x * squarel(rho[IDX(x, y+1, N)] - rho[IDX(x,y, N )]) + eq_activity[IDX(x ,y, N )].y * squarel(rho[IDX(x +1 , y,N )] - rho[IDX(x,y,N)]); 
		}
	}
	return diss/2; 
}

// The distance (discontinuous at the boundary )
//double distance( int x, int y, int a, int b ){
//	return ((x  -a)*(x- a  )+(y- b)*(y-b))/(N); 
//}
//
//// One Minimum
//void initialize_V(){
//	for(int i = 0; i < N; i++){
//		for( int j = 0; j < N ; j++){
//		V[IDX(i, j, N)] = distance(i,j , 30, 30 ) + distance(i, j,  70, 70 );
//		}				
//	}
//}
//
////// Distance (continuous at the boundary )
//double distance2 (int x, int y, int a, int b ){
//    	int dx = abs( x - a ); 
//	int dy = abs(y - b); 
//	dx = dx > N / 2 ? N - dx : dx;
//   	dy = dy > N / 2 ? N - dy : dy;
//	return dx + dy; 
//}
//
//void initialize2_V(){
//	for(int i = 0; i < N; i++){
//		for( int j = 0; j < N ; j++){
//		V[IDX(i, j, N)] = (distance(i, j, 30, 30 ))/100;
//		}				
//	}
//}
//void initialize3_V(){
//	for(int i = 0; i < N; i++){
//		for( int j = 0; j < N ; j++){
//		V[IDX(i, j, N)] = (abs(i - 10) *abs(j -10) + abs(i-50)*abs(j -60))/100;
//		}				
//	}
//}
//
//void initialize4_V() {
//
//	for(int i = 0; i < N; i++){
//		for( int j = 0; j < N ; j++){
//		V[IDX(i, j, N)] = fmin( distance2(i,j, 30, 30 ) , distance2(i, j, 70 , 70))/100;
//		}				
//	}
//}
//



//We take the rates (c^s(x,y) + c^a(x,y) ) e^{-\beta/2 V(y) - V(x)}, fo
// for an antisymmetric edge functiokn c^a- Usually c^s = 1, and c^s is small. \
// In order for e^{-\betaV}  to be the invariant measure, the following has to be satisfied 
// sum_{y \, y ~ x} c^a(x, y ) (e^{-\beta/2 V(y) } + e^{-\beta/2 V(x)}) = 0 
// The next function computes the above quantity 
// In the case of "orthogonal forces" this is clear, sicne c has itself divergence 0 and when $c$ is different from zero e^{-\beta/2  V(y)} + e^{-\beta/2 V} is constant.


// Here we mimich the case in which $V$ and $c$ are radial and orthogonal, and $c$ is of divergence 0
//void initialize5_V(){
//int max = (N +1) / 2;
//if (N % 2 != 1) {
//        fprintf(stderr, "Error: N must be odd, but N = %d\n", N);
//        return;  // exit with error code
//    }	
//// Define the N/2-1 different values of $V$
//double *a = malloc(((N+1)/2)* sizeof(double) ); 
//// Define the N/2 - 1 possible values fo $c$
//double *b = malloc((N+1)/2  * sizeof(double)); 
//	for(int i = 0; i < max; i++){
//	a[i ] = (double)(i*i)/N; // It took a lot of time to understand that the casting went wrong  
//	b[i] = sin (i); 
//	}
//	for(int i = 0; i < N* N; i++){
//	c_x[i] = 0; 
//	c_y[i ] = 0; 
//	}
//	V[IDX((N-1)/2 , (N - 1, N)/2) ] = a[0]; 
//	for(int i = 1; i <  max; i++ ){
//		int max2 = (N -1 )/2 +i ; 		
//		for(int j = (N-1)/2 - i  ; j < max2; j++){
//
//			V[IDX( (N-1)/2 + i , j , N)] = a[i];//Right side of the square  minus hte hifhrst point 
//			V[IDX((N-1)/2 - i, N - j- 1 , N)  ] = a[i];  	// left side of the square minus the lowest point 
//			V[IDX( N-1 - j, (N-1, N)/2 + i )] = a[i]; // upper sidei. The index that is boving is j 
//			V[ IDX(  j, (N-1, N)/2 - i) ] = a[i]; // lower side 
//		// c_x denotes the magnitude of the vector in the east dirsction, that is c((a,b ) , (a, b) + e_), while c_y in the north 
//			c_x[IDX( (N-1)/2 + i , j , N)] = -b[i]; //upper  side of the square 
//			c_x[IDX((N-1)/2 - i, j, N)  ] = b[i];  	// lower  side of the square
//			c_y[IDX( j , (N-1, N)/2 + i )] =  b[i]; // irght side
//			c_y[ IDX( j, (N-1, N)/2 - i) ] = -b[i];  // left side 
//			
//
//			}
//		
//	}
////			for (int i= N-1; i > -1; i--){
////				for (int j = 0; j <N ; j++ ){
////				printf("  (%f, %f)  " , c_x[IDX(i, j , N)], c_y[IDX(i, j, N)]); 
////				}
////				printf("\n"); 
////			}
////			for (int i= N-1; i > -1; i--){
////				for (int j =0; j < N; j++ ){
////				printf(" %f ", V[IDX(i, j, N)]); 
////				}
////				printf("\n"); 
////			}
//}



//void compute_mu() {
//    double Vmin = V[0];
//    for (int i = 1; i < N * N; i++) {
//        if (V[i] < Vmin) Vmin = V[i];
//    }
//
//    for (int i = 0; i < N * N; i++) {
//        mu[i] = exp(-beta * (V[i] - Vmin));  // Numerically stable
//    }
//
//}

// Drift field (example: zero)
//void drift(int x, int y, double *bx, double *by) {
//    *bx = 0.0;
//    *by = 0.0;
//}
//
//
//void normalize_mu() {
//    double Z = 0.0;
//    for (int i = 0; i < N * N; i++) Z += mu[i];
//    if (Z == 0.0) {
//        fprintf(stderr, "ERROR: Normalization factor Z is zero. Possibly overflow/underflow in mu.\n");
//        exit(1);
//    }
//    for (int i = 0; i < N * N; i++) mu[i] /= Z;
//}
//
//void compute_fields() {
//    compute_mu();
//    for (int y = 0; y < N; y++) {
//        for (int x = 0; x < N; x++) {
//            int xp = (x + 1) % N, xm = (x - 1 + N) % N;
//            int yp = (y + 1) % N, ym = (y - 1 + N) % N;
//
//            double dx = 0.5 * (V[IDX(xp, y, N)] - V[IDX(xm, y, N)]);
//            double dy = 0.5 * (V[IDX(x, yp, N)] - V[IDX(x, ym, N)]);
//            gradVx[IDX(x, y, N)] = dx;
//            gradVy[IDX(x, y, N)] = dy;
//
//            double bx, by;
//            drift(x, y, &bx, &by);
//            cx[IDX(x, y, N)] = bx + dx;
//            cy[IDX(x, y, N)] = by + dy;
//
//            mu[IDX(x, y, N)] = exp(-beta * V[IDX(x, y, N)]);
//        }
//    }
//    normalize_mu();
//    for (int i = 0; i < N * N; i++) {
//        mucx[i] = mu[i] * cx[i];
//        mucy[i] = mu[i] * cy[i];
//        divergence[i] = gradVx[i] * cx[i] + gradVy[i] * cy[i];
//    }
//}
void save_and_plot_vectors(Vector* vecs, int N, const char* name, double sigma, int skip, const char* python_script) {
    const char* binfile = "vector_frame.bin";

    // Save entire Vector array in one go
    FILE* f = fopen(binfile, "wb");
    if (!f) { perror("fopen"); exit(1); }
    fwrite(vecs, sizeof(Vector), N*N, f);
    fclose(f);

    // Build system command to call Python script
    char cmd[1024];
	snprintf(cmd, sizeof(cmd),
        	 "python3 %s %s %d %f %d '%s'",
        	 python_script, binfile, N, sigma, skip, name);    

    int ret = system(cmd);
    if (ret != 0) fprintf(stderr, "Python plotting failed!\n");
}



void save_and_plot_forms(Form* forms, int N, const char* label, int skip, const char* python_script) {
    // Create a unique temporary file
    char tmp_template[] = "/tmp/formsXXXXXX";
    int fd = mkstemp(tmp_template);
    if (fd == -1) {
        perror("Failed to create temporary file");
        return;
    }

    // Write the data
    FILE* f = fdopen(fd, "wb");
    if (!f) {
        perror("Failed to open temporary file stream");
        close(fd);
        return;
    }

    if (fwrite(forms, sizeof(Form), N * N, f) != (size_t)(N * N)) {
        perror("Failed to write data to file");
        fclose(f);
        return;
    }
    fclose(f);  // also closes fd

    // Prepare and execute the Python command
    char cmd[1024];
    snprintf(cmd, sizeof(cmd), "python3 %s %s %d %d \"%s\"", python_script, tmp_template, N, skip, label);
    int ret = system(cmd);
    if (ret != 0) {
        fprintf(stderr, "Python plotting failed!\n");
    }

    // Remove temporary file
    remove(tmp_template);
}


void plot_density_and_current(const char *dirname, int N_particles, int N) {
    // Build command safely
    char command[PATH_MAX * 2];
    snprintf(command, sizeof(command), 
             "python3 heatmap_empirical_measure_and_current.py --size %d --sigma 2 --particles %d "
             "--potential_type double_well --file \"%s/stats.bin\"",
             N, N_particles, dirname);

    // Execute command
    printf("Executing: %s\n", command);
    int ret = system(command);
    if (ret != 0) {
        fprintf(stderr, "Failed to run plot_current_and_density.py (exit code: %d)\n", ret);
        fprintf(stderr, "Command: %s\n", command);
    }
}

void plot_3d_density(const char *dirname, int N_particles, int N) {
    // Build command safely
    char command[PATH_MAX * 2];
    snprintf(command, sizeof(command), 
             "python3 animate_3d_density.py --file \"%s/stats.bin\" --particles %d "
             "--size %d --output \"%s\"",
             dirname, N_particles, N, dirname);

    // Execute command
    printf("Executing: %s\n", command);
    int ret = system(command);
    if (ret != 0) {
        fprintf(stderr, "Failed to run animate_3d_density.py (exit code: %d)\n", ret);
        fprintf(stderr, "Command: %s\n", command);
    }
}

void animate_particles(const char *dirname, int num_particles, int grid_size) {
    // Build command safely
    char command[PATH_MAX * 2];
    snprintf(command, sizeof(command), 
             "python3 animate_particles.py --file \"%s/trajectory.bin\" "
             "--particles %d --size %d --output \"%s\"",
             dirname, num_particles, grid_size, dirname);

    // Execute command
    printf("Executing: %s\n", command);
    int ret = system(command);
    if (ret != 0) {
        fprintf(stderr, "Failed to run animate_particles.py (exit code: %d)\n", ret);
        fprintf(stderr, "Command: %s\n", command);
    }
}
