// init.h

#ifndef INIT_H
#define INIT_H

extern int global_variable;


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
