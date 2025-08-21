#ifndef INPUT_OUTPUT_H
#define INPUT_OUTPUT_H

#include <stddef.h>  // for size_t

#define MAX_LINE 256

/// Configuration contains the input required by the program
typedef struct {
    int N;              // scaling parameter
    double T;           // Macroscopic time 
    double L;           // Macroscopic length
    int N_particles;    // Number of particles
    double beta;        // Inverse temperature 
    int resolution;     // The discretization step 
} Configuration;

typedef int (*ConfigSetter)(Configuration *conf, const char *value);

typedef struct {
    const char *key;
    ConfigSetter setter;
} ConfigEntry;

// Setter functions
int set_N(Configuration *conf, const char *value);
int set_T(Configuration *conf, const char *value);
int set_L(Configuration *conf, const char *value);
int set_N_particles(Configuration *conf, const char *value);
int set_resolution(Configuration *conf, const char *value);
int set_beta(Configuration *conf, const char *value);

// Core config functions
void read_config_file(Configuration *conf, const char *filename);
int is_comment_or_blank(const char *line);

// Output helpers
void print_config_to_screen(Configuration conf);
void create_output_directory(const Configuration conf, char *dirname_out, size_t len);
void save_config_to_json(const Configuration conf, const char *dirname);

#endif // INPUT_OUTPUT_H

