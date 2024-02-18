//terminate_block.h

#ifndef TERMINATE_BLOCK_H
#define TERMINATE_BLOCK_H


/*------------------GLOBAL VARIABLES-----------------------------------------*/
/*---------------------------------------------------------------------------*/

/*------------------STRUCTURE DEFINITIONS------------------------------------*/
/*---------------------------------------------------------------------------*/

/*------------------GLOBAL STRUCTURES----------------------------------------*/
/*---------------------------------------------------------------------------*/

/*------------------GLOBAL FUNCTIONS-----------------------------------------*/
extern void print_input(Slice *);
extern void conf_output(Slice *);
extern void printenergy(Slice *);
extern void printenergy_warmup(Slice *);
extern void print_slice_information(Slice *);
extern void printing_trajectory(Slice *);
extern void append_to_file(char [100], int , double );
extern void print_to_file(char [100], int , double );
extern void append_to_file2(char [100], char [1000] , char [100]  );
extern void print_to_file2(char [100], char [1000] , char [100] );
extern void terminate_block(Slice *);
extern void finalstat(Slice *);
/*---------------------------------------------------------------------------*/

#endif
