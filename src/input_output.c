// Here we write the functions used to read the input parameters and to write the outputs 



#include <stdio.h>      // snprintf, printf, etc.
#include <stdlib.h>     // exit, malloc, free, etc.
#include <string.h>     // strcmp, strcpy, etc.
#include <ctype.h>      // isspace
#include <sys/stat.h>   // mkdir
#include <sys/types.h>  // mode_t (needed on some systems)
#include <errno.h>      // errno, EEXIST
#include "input_output.h"


ConfigEntry config_table[] = {
    {"N", set_N},
    {"T", set_T},
    {"L", set_L},
    {"N_particles", set_N_particles}, 
    {"resolution", set_resolution},
    {"beta", set_beta}, 
    {NULL, NULL} //A sentinel indicating the end of the structure. 
}; 
// This associates to the "key", the first element,  a function set_config, the second entry.  
// Similarly, in the .txt file to read, to each key there is associated a value, which we would like to be assigned to the corresponding variable in the program.  
// To this end, we need to choose a precise function. This table does the job: the function associated to the key assigns values to the corresponding variable in teh program  



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////  OUTPUT DIRECTORIES AND PLOTS ///////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////






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




void print_config_to_screen(Configuration conf){	
 printf("You have made the following choices: \n");
 printf("   N: %d,\n\n",conf.N);
 printf("   T: %f \n\n", conf.T);
 printf("   L: %f \n\n", conf.L);
 printf("   N_particles: %d \n\n", conf.N_particles);
 printf("   Resolution: %d,\n\n", conf.resolution);
 printf("   Beta =  %f,\n\n", conf.beta); 
}

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


