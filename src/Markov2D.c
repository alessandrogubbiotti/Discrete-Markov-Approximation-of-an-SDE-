////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// 2025-08-15 Main changes that need to be done: 
//// 0) A good normalization for the empirical current, add beta to the output directory name
//// 1) Write the parameters in the json output file and make the scripts read from there. Then ask to the user 
//// if he wishes to save the results or visualize them 
//// 2) Add an estimation of c and V and compare it with norm infty, 1 and 2, with our actual fields. 
//// 3) Idem for the invariant measure. This is in light to do it when we do not know the quantities 
//// 4) If possible, estimate the mixing time and see that the non-reversible dynamics is faster
////
//// What I am going to do next is to estimate the dissipation potentials 
//// In a totally different perspective, I should also write a stream field, maybe using the other program that
//// I have made. 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// THE INPUT VARIABLES ARE 
// N The number of sites for each dimension. We consider an NxN grid with periodic boundary conditions
// L The macroscopic length of the grid. It has been not used  untill now 
// T Macroscopic time
// N_particles is the number of particles in the simulation
// 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// THE RATES ARE OF THE FORM 
// (1 + c(a,b)) e^{ -beta/2 (V(b) - V(a)}, 
// for each a and b in the grid. and c antisymmetric 
// and we assume that sum_{y} c(x,y) e^{ -\beta/2(V(y) + V(x))}  = 0
// Untill now, the only set of examples that I have are of the form 
// c is divergence free and orthogonal to V in the following sense  
// for each a, b in the grid, either V(b) - V(a) = 0 or c(a,b) = 0
//////////////////////////////////////////////////////////////////////////
// To save  funcitons of the directed edges we use the struct Vector
// To sace antisymmetric functions of the edges we use Form and we take the convention to 
// that along the x axis we take the sign that the funciton has on the edge in the east direction, and for the
// y-axis along the north direction
// To obtain the south direction, for instance, we should look at the point one unit to the south and take the 
// negative of its y-component.
// In this way we are able to plot currents as vectors 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// Since a diffusive rescaling should lead to an Ito diffusion for each single particle
//// we take a macroscopic resolution in which the macroscopic unit of time,
//// composed by N*N units of micrscopic time, in resolution subintervals (Have to change this now)
//// Therefore
////
////
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
#include "analysis.h"
#include "Markov2D.h"
#include "input_output.h"

#define MAX_LINE 256




gsl_rng *rng; ///GSLsetyup




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// HEADERS ////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void current_distances(Form *a, Form *b, int N);
void typical_current(Form *J_typ, Vector *rates, double *empirical_measure,  double time_step, int N  ); 

void set_skew_symm_current(Form *skew_symm_current, long double *rho, Form *eq_current,  double dt, int N); 

void diagnose_kinetic_cost_variants(Form *empirical_current,
                                    Form *eq_current,
                                    Activity *eq_activity,
                                    long double *rho,
                        		double dt,   
			          int N);


void normalize_empirical_flow(Vector *empirical_flow, int resolution,  int N_particles,  int NN  );
void normalize_empirical_measure(double *empirical_measure, int resolution,  int N_particles,  int NN  ); 
void diagnose_currents(long double *m_beta, long double *rho,  Form *empirical_current, Form *eq_current, int N ); 

void divide_vector(Vector *vec, double scalar, int N );
void divide_form(Form *curr, long double scalar, int N);
 
void plot_3d_density(const char *dirname, int N_particles, int N);
void animate_particles(const char *dirname, int num_particles, int grid_size);
void plot_density_and_current(const char *dirname, int N_particles, int N);
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// INPUT READING  /////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////



 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////// READING THE INPUTS /////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



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
long double kinetic_cost(Form *j, Activity *eq_activity, Form *eq_current, long double *rho, double time,  int  N  ); 
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
//	save_and_plot_vectors(rates_flat, N, "rates", 0.05, 2, "plot_vectors.py" ); 
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
	Form *skew_symm_current = malloc(N*N*sizeof(Form)); 
	Form *null_current = calloc(N*N, sizeof(Form)); 
	invariant_measure( m_beta, V, beta , N); 
	equilibrium_activity(eq_activity, m_beta, rates_flat, N); 
	equilibrium_current(eq_current, m_beta, rates_flat, N); 
   	
//	save_and_plot_forms( dV, N, "GradV", 2, "plot_forms.py"); 
//	save_and_plot_forms( c, N, "c", 2, "plot_forms.py"); 


	
	////////////////////////////////////////////////////////
	///////////// SIMULATION AND STATISTICS ////////////////
	////////////////////////////////////////////////////////

int count = 0; 

	
	// We start by allocating the empirical quantities we will use 
	//DEBUG TYPICAL CURRENT 
	Form *J_typ = malloc(N*N*sizeof(Form)); //typical_current 
	Form *empirical_current = malloc(N*N*sizeof(Form)); 
	Vector *empirical_flow = malloc(N*N*sizeof(Vector)); 
	Vector *cumulative_empirical_flow = malloc(N*N*sizeof(Vector)); 
	double *empirical_measure = malloc(N*N*sizeof(double)); 	
	double *cumulative_empirical_measure = malloc(N*N*sizeof(double)); 	
	
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
		
// DIAGNOSTIC: compare K(a) for q = eq_current * (rho_i + rho_j) and K_true for J(mu)
long double S_j2 = 0.0L, S_jq = 0.0L, S_q2 = 0.0L;
long double comp_j2 = 0.0L, comp_jq = 0.0L, comp_q2 = 0.0L;
long double S_j2_true = 0.0L, S_jq_true = 0.0L, S_q2_true = 0.0L; // not used but kept
long double total_j2 = 0.0L, comp_total = 0.0L;
int a = (N*N)/resolution; 
for (int x = 0; x < N; ++x) {
  for (int y = 0; y < N; ++y) {
    int idx = IDX(x,y,N);
    // horizontal edge: (x,y) -> (x, y+1)
    int idx_r = IDX(x,(y+1)%N,N);
    long double mu_i = (long double) empirical_measure[idx];
    long double mu_j = (long double) empirical_measure[idx_r];

    // exact expected current on that edge (J(mu)):
    long double rate_east_i = (long double) rates[x][y].east;
    long double rate_west_j = (long double) rates[x][(y+1)%N].west;
    long double Jmu_x = (mu_i*rate_east_i - mu_j*rate_west_j)*a;

    // your approximate q: eq_current * (rho_i + rho_j)
    long double qx = eq_current[idx].x * (rho[idx] + rho[idx_r])*a; // if you multiply by alpha later

    long double jx = empirical_current[idx].x;
    long double Ax = eq_activity[idx].x > 0.0L ? eq_activity[idx].x : 1e-300L;

    kahan_add(jx*jx/Ax, &S_j2, &comp_j2);
    kahan_add(jx*qx/Ax, &S_jq, &comp_jq);
    kahan_add(qx*qx/Ax, &S_q2, &comp_q2);

    // also accumulate true residual sums if you want to inspect (with Jmu)
    kahan_add(jx*jx/Ax, &S_j2_true, &comp_total);
    kahan_add(jx*Jmu_x/Ax, &S_jq_true, &comp_total);
    kahan_add(Jmu_x*Jmu_x/Ax, &S_q2_true, &comp_total);

    // vertical edge: (x,y) -> (x+1,y)
    int idx_u = IDX((x+1)%N,y,N);
    mu_i = (long double) empirical_measure[idx];
    mu_j = (long double) empirical_measure[idx_u];
    long double rate_north_i = (long double) rates[x][y].north;
    long double rate_south_j = (long double) rates[(x+1)%N][y].south;
    long double Jmu_y = (mu_i*rate_north_i - mu_j*rate_south_j)*a;

    long double qy = eq_current[idx].y * (rho[idx] + rho[idx_u])*a;
    long double jy = empirical_current[idx].y;
    long double Ay = eq_activity[idx].y > 0.0L ? eq_activity[idx].y : 1e-300L;

    kahan_add(jy*jy/Ay, &S_j2, &comp_j2);
    kahan_add(jy*qy/Ay, &S_jq, &comp_jq);
    kahan_add(qy*qy/Ay, &S_q2, &comp_q2);

    kahan_add(jy*jy/Ay, &S_j2_true, &comp_total);
    kahan_add(jy*Jmu_y/Ay, &S_jq_true, &comp_total);
    kahan_add(Jmu_y*Jmu_y/Ay, &S_q2_true, &comp_total);
  }
}
// algebra for q-ansatz
long double a_star = (S_q2 > 0.0L) ? (S_jq / S_q2) : 0.0L;
printf("DIAG_S: S_j2=%.12Lg  S_jq=%.12Lg  S_q2=%.12Lg  a* = %.6Lg\n",
       S_j2, S_jq, S_q2, a_star);

// Print predicted K for alphas you used:
long double alphas[] = {0.0L, 0.25L, 0.5L, 1.0L};
for (int ia = 0; ia < 4; ++ia) {
  long double a = alphas[ia];
  long double S_a = S_j2 - 2.0L*a*S_jq + a*a*S_q2;
  long double K_a = 0.5L * S_a;
  printf("predicted K(alpha=%.2Lf) = %.12Lg\n", a, K_a);
}

// Also compute K_true where we subtract exact J(mu):
long double S_true = S_j2_true - 2.0L*S_jq_true + S_q2_true;
long double K_true = 0.5L * S_true;
printf("K_true (subtract exact J(mu)) = %.12Lg\n", K_true);


		diagnose_kinetic_cost_variants( empirical_current, eq_current, eq_activity, rho, (double) resolution,  N); 

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
		
		normalize_empirical_measure(empirical_measure, resolution, N_particles, N*N ); 
		density(rho, m_beta, empirical_measure, N );

		normalize_empirical_flow(empirical_flow, N_particles, resolution,  N*N); 
		flow_to_current(empirical_current, empirical_flow,N); 	
		

		

		entropy = relative_entropy( rho, m_beta, N); 
		cost += kinetic_cost( empirical_current, eq_activity, eq_current, rho , N*N/resolution , N) * (N*N)/((long double)resolution); 
		potential += dissipation_potential( eq_activity, rho, N) * (N*N)/((long double)resolution); 
		fprintf(dissipation, "%f %.5Lf %.5Lf %.5Lf\n",  t_macro, entropy, cost, potential); 		
		
		set_skew_symm_current(skew_symm_current, rho,eq_current, (N*N)/(double)resolution ,N); 
		typical_current( J_typ,rates_flat,  empirical_measure, /*1/((double) resolution)*/ N*N /(double)resolution,  N  ); 
		current_distances(empirical_current, skew_symm_current, N );
		if( count == 10 ){
			count = 0;  
			save_and_plot_forms(empirical_current, N, "empirical current", 2 , "plot_forms.py" ); 
			save_and_plot_forms(J_typ, N, "typical current  current", 2 , "plot_forms.py" ); 
			save_and_plot_forms(skew_symm_current, N, "skew symmetric current", 2 , "plot_forms.py" ); 
			save_and_plot_forms(eq_current, N, "equilibrium current", 2 , "plot_forms.py" ); 
			double rho_sum = 0; 
			for(int a = 0; a < N; a++ ){
				for(int b = 0; b < N; b++){
					rho_sum += m_beta[IDX(a,b,N)] * rho[IDX(a,b, N )]; 
				}
			}
			printf( "\n\n The total sum of m_beta * rho is  %f \n\n", rho_sum); 
			current_distances(empirical_current, null_current, N ); 
		}
		count++; 
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
	free(J_typ); free(skew_symm_current); free(null_current);
	free(cumulative_empirical_flow); free(cumulative_empirical_measure);   
	/////////////////////////////////////////////////////////////////////////
	///////////// SYSTEM CALLS TO PLOT STUFF ////////////////////////////////
	/////////////////////////////////////////////////////////////////////////
	
//	plot_density_and_current(dirname, N_particles, N);
//       plot_3d_density(dirname, N_particles, N); 
//	animate_particles(dirname, N_particles , N);



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
	int t_simulation = N*N/resolution;  
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
		  	if(  gsl_rng_uniform(rng) < 1 -  exp( - rate /(double)(resolution) ) ){ // We only move when the clock rings 
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
 	      if( i < N_particles/2 ){
	       particles[i].x = (N-1)/2;
 	       particles[i].y = (N-1)/4;
		}
		else {
	       particles[i].x = (N-1)/2;
 	       particles[i].y = 3*(N-1)/4;		
		}
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




//////////////////////////////////////////////////////////////////
/////// FUNCTIONS FOR READING THE INPUT PARAMETERS ///////////////
//////////////////////////////////////////////////////////////////




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
    const int NN = N * N;
    double Vmin = 0.0;
    long double S = stable_partition_sum(V, beta, NN, &Vmin); // Kahan + shift

    // m_beta[i] = exp(-beta (V[i] - Vmin)) / S  -> sums to 1
    for (int i = 0; i < NN; ++i) {
        m_beta[i] = expl(-(long double)beta * ((long double)V[i] - (long double)Vmin)) / S;
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
			eq_curr[IDX(x,y,N)].y = (m_beta[IDX(x,y,N)] * rates[IDX(x,y, N)].north -  m_beta[IDX(x+1,y,N)] * rates[IDX(x +1,y,N)].south); 
			eq_curr[IDX(x,y,N)].x = (m_beta[IDX(x,y,N)] * rates[IDX(x,y,N)].east - m_beta[IDX(x,y+1,N)] * rates[IDX(x,y+1, N)].west); 
		}
	}
}
void equilibrium_activity(Activity *eq_activity, long double *m_beta, Vector *rates, int N  ){
	for(int x = 0; x < N; x++){
		for(int y = 0; y < N; y++){
			eq_activity[IDX(x, y, N )].y = (m_beta[IDX(x,y,N)] * rates[IDX(x,y, N)].north  +  m_beta[IDX(x+1,y,N)] * rates[IDX(x +1,y,N)].south)/2.0L;
			eq_activity[IDX(x,y,N)].x = (m_beta[IDX(x,y,N)] * rates[IDX(x,y,N)].east  +  m_beta[IDX(x,y+1,N)] * rates[IDX(x,y+1, N)].west)/2.0L;
		}
	}
}




////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
//////////////  FLAT STRUCTURES  ON THE TANGENT SPACES  ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////


void  d(Form *dV, double *V, int N){
	for(int x = 0; x < N; x++){
		for(int y = 0; y < N; y++ ){
			dV[IDX(x,y, N)].x = V[IDX(x  , y +1, N )]  - V[IDX(x,y,N)]; 
			dV[IDX(x,y,N)].y = V[IDX(x +1,y, N )] - V[ IDX(x,y,N)]; 
		}
	} 
}



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
	return (r - 1)*(r - 1); 
}

/// Il dt sta fuori o sta dentro??????
long double kinetic_cost(Form *j, Activity *eq_activity, Form *eq_current, long double *rho,  double  dt, int N){
    long double sum = 0.0L, comp = 0.0L;

    for(int x = 0; x < N; x++){
        for(int y = 0; y < N; y++){
            // horizontal edge (x, y) -> (x, y+1)
            long double avg_rho_x = 0.5L *  (rho[IDX(x,y,N)] + rho[IDX(x, y+1, N)]);
            long double res_x = j[IDX(x,y,N)].x - eq_current[IDX(x,y,N)].x * avg_rho_x * dt;
            long double Ax = eq_activity[IDX(x,y,N)].x;
            long double term_x = (Ax > 0.0L) ? (res_x*res_x)/Ax : 0.0L;

            // vertical edge (x, y) -> (x+1, y)
            long double avg_rho_y = 0.5L * (rho[IDX(x,y,N)] + rho[IDX(x+1, y, N)]);
            long double res_y = j[IDX(x,y,N)].y - eq_current[IDX(x,y,N)].y * avg_rho_y * dt ;
            long double Ay = eq_activity[IDX(x,y,N)].y;
            long double term_y = (Ay > 0.0L) ? (res_y*res_y)/Ay : 0.0L;

            kahan_add(term_x + term_y, &sum, &comp);
        }
    }
    return 0.5L * sum;
}

// Dissipation potential consistent with the inequality for quadratic entropy
// where h(r) = (r-1)^2, so h'(r) = 2(r-1) and Delta h' = 2 * (r_j - r_i).
// The dissipation per edge is (1/2) * A_eq * |Delta h'|^2 = 2 * A_eq * (r_j - r_i)^2.
//
// Uses Kahan compensated summation to reduce rounding error.
long double dissipation_potential(Activity *eq_activity, long double *rho, int N){
    long double sum = 0.0L, comp = 0.0L;
    for (int x = 0; x < N; ++x) {
        for (int y = 0; y < N; ++y) {
            int idx = IDX(x,y,N);

            // east edge: (x,y) -> (x, y+1)
            long double Ax = eq_activity[idx].x;
            long double r_i = rho[idx];
            long double r_j = rho[IDX(x, y+1, N)];
            long double dr = r_j - r_i;
            long double term_x = 0.0L;
            if (Ax > 0.0L) {
                // (1/2) * A * (Delta h')^2 = (1/2) * A * (2 dr)^2 = 2 * A * dr^2
                term_x = 2.0L * Ax * dr * dr;
            }

            // north edge: (x,y) -> (x+1, y)
            long double Ay = eq_activity[idx].y;
            long double r_i2 = rho[idx];
            long double r_j2 = rho[IDX(x+1, y, N)];
            long double dry = r_j2 - r_i2;
            long double term_y = 0.0L;
            if (Ay > 0.0L) {
                term_y = 2.0L * Ay * dry * dry;
            }

            kahan_add(term_x + term_y, &sum, &comp);
        }
    }
    return 0.5L*sum; // note: already summed with the 1/2*(Delta h')^2 factor -> returns full D
}


//
//long double dissipation_potential(Activity *eq_activity, long double *rho, int N){
//    long double sum = 0.0L, comp = 0.0L;
//
//    for(int x = 0; x < N; x++){
//        for(int y = 0; y < N; y++){
//            long double Ax = eq_activity[IDX(x,y,N)].x;
//            long double Ay = eq_activity[IDX(x,y,N)].y;
//
//            long double drx = rho[IDX(x, y+1, N)] - rho[IDX(x, y, N)];
//            long double dry = rho[IDX(x+1, y, N)] - rho[IDX(x, y, N)];
//
//            long double term = (Ax > 0.0L ? Ax * drx * drx : 0.0L)
//                             + (Ay > 0.0L ? Ay * dry * dry : 0.0L);
//
//            kahan_add(term, &sum, &comp);
//        }
//    }
//    return 0.5L * sum;
//}
//

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




void divide_vector(Vector *vec, double scalar, int N ){
	for(int i= 0; i < N*N; i++){
		vec[i].north = vec[i].north/scalar; 
		vec[i].south = vec[i].south/scalar; 
		vec[i].east = vec[i].east/scalar; 
		vec[i].west = vec[i].west/scalar; 
	}
}

void divide_form(Form *curr, long double scalar, int N){
	for(int i = 0; i < N*N; i++){
		curr[i].x = curr[i].x/scalar; 
		curr[i].y = curr[i].y/scalar; 
	}
}

//
//
//void plot_entropy_dissipation_inequlaity(const char *dirname) {
//    // Build command safely
//    char command[PATH_MAX * 2];
//    snprintf(command, sizeof(command), 
//             "python3 dissipation.py");
//
//    // Execute command
//    printf("Executing: %s\n", command);
//    int ret = system(command);
//    if (ret != 0) {
//        fprintf(stderr, "Failed to run animate_particles.py (exit code: %d)\n", ret);
//        fprintf(stderr, "Command: %s\n", command);
//    }
//}



//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// QUANTITIELS FUNCTIONS THAT QUANTIFY THE CONVERGENCE/////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// DEBUG  /////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

// assume: NN = N*N, total_time (double), c (long double), m_beta[i] (long double)
// arrays: Nplus_x[idx], Nminus_x[idx] for +x and -x transitions (or general edges)
// for simplicity I'll show horizontal edges only; repeat for vertical.

void diagnose_currents(long double *m_beta, long double *rho,  Form *empirical_current, Form *eq_current, int N ) {
	long double  NN = N*N;
	long double sum_weighted_res2 = 0.0L, comp = 0.0L; // Kahan
	long double max_abs_res = 0.0L;
	int max_idx1 = -1;
	int max_idx2 = -1; 
	for (int x = 0; x < N; x++) {
		for( int y = 0; y < N; y++){	 
			long double mi = m_beta[IDX(x,y,N)]/NN;
			long double mj = m_beta[IDX(x, y+1, N )]/NN;
			long double s = +1.0L; // for +x edge
			//long double J_eq  = 2.0L * s * c * sqrtl(mi*mj);  // must match your m_beta scale
			long double res = empirical_current[IDX(x,y,N)].x - 0.5L*eq_current[IDX(x,y,N)].x*(rho[IDX(x,y,N)] + rho[IDX(x,y+1,N)]);
			
			// long double Aeq = /* your activity value for that edge, must be >0 */ sqrtl(mi*mj); // example
			// weighted residual^2 / Aeq
			long double Aeq = 2*sqrtl(mi*mj); 
			long double term = (res*res) / ( Aeq > 0.0L ? Aeq : 1e-300L );
			// Kahan add
			long double k = term - comp;
			long double t = sum_weighted_res2 + k;
			comp = (t - sum_weighted_res2) - k;
			sum_weighted_res2 = t;

      			if (fabsl(res) > max_abs_res) { 
				max_abs_res = fabsl(res); 
				max_idx1 = x; 
				max_idx2 = y;
			 }
    		}
	}
	long double total_cost = sum_weighted_res2;
	printf("diagnose: total_kinetic_cost (un-normalized) = %.12Lg\n", total_cost);
	printf("diagnose: max_abs_res = %.6Lg at idx (%d, %d)\n", max_abs_res, max_idx1, max_idx2);
}



void current_distances(Form *a, Form *b, int N) {
    int NN = N*N;
    long double l1 = 0.0L, l2sq = 0.0L, linf = 0.0L;
    for (int i = 0; i < NN; i++) {
        long double dx = a[i].x - b[i].x;
        long double dy = a[i].y - b[i].y;
        long double norm = sqrtl(dx*dx + dy*dy); // Euclidean norm at each site
        l1   += fabsl(norm);
        l2sq += norm*norm;
        if (norm > linf) linf = norm;
    }
    printf("The L^2 distance is: %Lf\n  The L^infty distance is %Lf\n	The L^1 distance is %Lf\n\n",l2sq, linf,  l1); 
}


void typical_current(Form *J_typ, Vector *rates, double *empirical_measure, double time_step,  int N  ){
	for(int x = 0; x < N; x++){
		for(int y = 0; y < N; y++){
			J_typ[IDX(x,y,N)].x =  (empirical_measure[IDX(x,y,N)] * rates[IDX(x,y,N)].east - empirical_measure[IDX(x,y+1, N)] * rates[IDX(x,y+1,N)].west)*time_step;
			J_typ[IDX(x,y,N)].y =  (empirical_measure[IDX(x,y,N)] * rates[IDX(x,y,N)].north - empirical_measure[IDX(x+1,y, N)] * rates[IDX(x+1,y,N)].south)*time_step;
		}
	}
}


void normalize_empirical_measure(double *empirical_measure, int resolution,  int N_particles,  int NN  ){
	int a = NN/resolution; 
	for (int i = 0; i < NN; i++){
		empirical_measure[i] = (empirical_measure[i])/((double)(N_particles* a) ); 
	}
}

void normalize_empirical_flow(Vector *empirical_flow,  int N_particles,  int dt, int NN  ){
	for (int i = 0; i < NN; i++){
		empirical_flow[i].north = (empirical_flow[i].north)/((double)N_particles); 
		empirical_flow[i].south = (empirical_flow[i].south)/((double)N_particles); 
		empirical_flow[i].east = (empirical_flow[i].east)/((double)N_particles); 
		empirical_flow[i].west = (empirical_flow[i].west)/((double)N_particles ); 
	}
}

long double relative_entropy(long double *rho, long double *m_beta, int N){
    long double sum = 0.0L, comp = 0.0L;

    for (int i = 0; i < N*N; i++){
        long double term = entropy(rho[i]) * m_beta[i];
        kahan_add(term, &sum, &comp);
    }

    return sum;
}


void diagnose_kinetic_cost_variants(Form *empirical_current,
                                    Form *eq_current,
                                    Activity *eq_activity,
                                    long double *rho,
					double time,  
                                    int N)
{
    long double alphas[] = {0.0L,0.25L, 0.5L, 1.0L};
    int n_alphas = 4;

    for (int a = 0; a < n_alphas; ++a) {
        long double alpha = alphas[a];
        long double total_cost = 0.0L;
        long double comp = 0.0L; // for Kahan summation

        for (int x = 0; x < N; ++x) {
            for (int y = 0; y < N; ++y) {
                int idx = IDX(x,y,N);

                // Horizontal edge (x,y) -> (x, y+1)
                int idx_right = IDX(x, (y+1)%N, N);
                long double qx = alpha * eq_current[idx].x *
                                 (rho[idx] + rho[idx_right])*((double)N*N) /((double)time) ;
                long double jx = empirical_current[idx].x;
                long double Ax = eq_activity[idx].x > 0.0L ? eq_activity[idx].x : 1e-300L;
                long double resx = jx - qx;
                kahan_add((resx*resx)/Ax, &total_cost, &comp);

                // Vertical edge (x,y) -> (x+1, y)
                int idx_up = IDX((x+1)%N, y, N);
                long double qy = alpha * eq_current[idx].y *
                                 (rho[idx] + rho[idx_up])*((double)N*N) /((double)time);
                long double jy = empirical_current[idx].y;
                long double Ay = eq_activity[idx].y > 0.0L ? eq_activity[idx].y : 1e-300L;
                long double resy = jy - qy;
                kahan_add((resy*resy)/Ay, &total_cost, &comp);
            }
        }

        total_cost *= 0.5L; // match your kinetic_cost normalization
        printf("Kinetic cost with alpha= - %.2Lf : %.12Lg\n",  alpha, total_cost);
	}
}
		

void set_skew_symm_current(Form *skew_symm_current, long double *rho, Form *eq_current,  double dt, int N){
	for(int x = 0; x < N; x++){
		for(int y = 0; y < N; y++){
			skew_symm_current[IDX(x,y,N )].x = 0.5L*(rho[IDX(x,y,N)] + rho[IDX(x,y+1, N)])*eq_current[IDX(x,y,N)].x*dt; 
			skew_symm_current[IDX(x,y,N )].y = 0.5L*(rho[IDX(x,y,N)] + rho[IDX(x + 1, y, N)])*eq_current[IDX(x,y,N)].y*dt ; 
		}
	}
}




