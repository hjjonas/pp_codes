// init.h

#ifndef INIT_H
#define INIT_H

// extern int global_variable;


/*---------------------------------------------------------------------------*/
/*------------------LOCAL FUNCTIONS------------------------------------*/

void check_input_with_MAXDEFS(void);
vector check_read_unitvec(vector );
double find_forcecutoff_distance(void );
double find_minimum_of_potential(void);
double find_trunc_of_Saccent_1(int);
double find_zcut(double );
double find_overlap_distance(void);
void gravitational_parameters(void);
void init_model(Slice *);
void init_simtype(Slice *);
void print_input(Slice *);
void print_particle_properties(void);
void ratcheting_prep(Slice *);
int ratcheting_prep_conditions(int , int );
void read_input(Slice *);
void restart_read_time(Slice *);
void setup_delta(void);
void setup_mc_move(MC *);


/*---------------------------------------------------------------------------*/
/*------------------STRUCTURE DEFINITIONS------------------------------------*/

typedef struct init_type {
// Here all parameters to initialize the configuration, e.g. read from file 
    
    int           mc_warmup,
                  empty_files,
                  nchains,
                  restart,
                  read_path,      //  this is used if you want to read trajectory.xyz/inp 
                  start_type,
                  particle_setup;

    double        chaingap;      // the distance between two chains when putting them in a row

    char          directorypath[500];   // used to load files in a different folder, else this is "./" i.e. the current folder
} Init;


#endif
