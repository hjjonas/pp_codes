//terminate_block.h

#ifndef TERMINATE_BLOCK_H
#define TERMINATE_BLOCK_H


/*------------------GLOBAL VARIABLES-----------------------------------------*/
/*---------------------------------------------------------------------------*/

#define MAX_FILENAME_LENGTH 100
#define MAX_EXT_LENGTH 100
#define MAX_VALUE_LENGTH 1000
/*------------------STRUCTURE DEFINITIONS------------------------------------*/
/*---------------------------------------------------------------------------*/

/*------------------GLOBAL STRUCTURES----------------------------------------*/
/*---------------------------------------------------------------------------*/

/*------------------GLOBAL FUNCTIONS-----------------------------------------*/
extern void print_input(Slice *);
extern void conf_output(Slice *);
extern void printenergy(Slice *, char [MAX_EXT_LENGTH]);
extern void print_slice_information(Slice *);
extern void print_StatsLength_to_file(StatsLength *);
extern void printing_trajectory(Slice *);
extern void write_append_to_file(char [100],  char [100], char ,char [1000] );
extern void terminate_block(Slice *);
extern void finalstat(Slice *);
/*---------------------------------------------------------------------------*/

#endif
