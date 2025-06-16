// version_specific_functions.h
#ifndef vsf_H
#define vsf_H

/*------------------GLOBAL VARIABLES-----------------------------------------*/
/*---------------------------------------------------------------------------*/
#define MAXBIN 200


/*---------------------------------------------------------------------------*/

/*------------------GLOBAL STRUCTURES-----------------------------------------*/
typedef struct LocalDensitytype  {
  // the 1D array of the local density  
  // the max density is approx 1 , so local_density \in [0,2 ]     

  StatsLength    local_density;
  
  double         ld_gridsize,
                 box_gridarea; // 

} LocalDensity ;

LocalDensity ld_array[MAXBIN];

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
