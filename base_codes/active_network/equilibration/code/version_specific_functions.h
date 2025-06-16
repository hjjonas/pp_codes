// version_specific_functions.h

#ifndef vsf_H
#define vsf_H

/*------------------GLOBAL VARIABLES-----------------------------------------*/
/*---------------------------------------------------------------------------*/
#define MAXBIN 100

/*------------------STRUCTURE DEFINITIONS------------------------------------*/
typedef struct histogram_type {
  /* historgam*/

    long        bin[NPART][MAXBIN], // first dimension is red bond number and other is bin of thistogram
                length[MAXBIN];

    double      min,
                max,
                d_bin;

    int         nbins;

} Histogram;


/*---------------------------------------------------------------------------*/

/*------------------GLOBAL STRUCTURES-----------------------------------------*/
extern Histogram  histogram_S, histogram_E, histogram_r;
extern int reduced_bondnumbers[NPART][NPART], Nred_bonds;
/*---------------------------------------------------------------------------*/


/*------------------GLOBAL FUNCTIONS-----------------------------------------*/
extern void innerloop_analysis(Slice *);
extern void version_specific_analysis(Slice *);
extern void every_timestep_analysis(Slice *);

extern void special_init(Slice *);
extern void allocate_memory(Slice *psl);
extern void free_all_memory(void);


/*---------------------------------------------------------------------------*/




#endif
