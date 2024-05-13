// version_specific_functions.h

#ifndef vsf_H
#define vsf_H

/*------------------GLOBAL VARIABLES-----------------------------------------*/
/*---------------------------------------------------------------------------*/

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
/*---------------------------------------------------------------------------*/


/*------------------GLOBAL FUNCTIONS-----------------------------------------*/
extern void innerloop_analysis(Slice *);
extern void version_specific_analysis(Slice *);
extern void every_timestep_analysis(Slice *);

extern void initialize_bond_configurations_measurements(Slice *);
/*---------------------------------------------------------------------------*/




#endif
