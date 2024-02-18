//terminate_block.h

#ifndef TERMINATE_BLOCK_H
#define TERMINATE_BLOCK_H

/*---------------------------------------------------------------------------*/
/*------------------GLOBAL VARIABLES-----------------------------------------*/


/*---------------------------------------------------------------------------*/
/*------------------STRUCTURE DEFINITIONS------------------------------------*/


/*---------------------------------------------------------------------------*/
/*------------------LOCAL FUNCTIONS------------------------------------------*/

void printstatusbmd(Slice *psl);
void count_empty_clusters_monomers(Slice *, int *, int *);
void print_pos_sites(Slice *);
void print_energy_s_theta1_theta2(Slice *);
void print_statistics_N1(Statistics , char [100]);
void print_statistics_file(StatsLength *,Slice *);
void printrdf();
// PUT IN VERSION SPECIFIC c file: void print_association_dissociation_times(void);

/*---------------------------------------------------------------------------*/
/*------------------GLOBAL STRUCTURES-----------------------------------------*/

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



#endif
