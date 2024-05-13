// bmd.h

#ifndef BMD_H
#define BMD_H

/*------------------STRUCTURE DEFINITIONS------------------------------------*/
/*---------------------------------------------------------------------------*/

/*------------------GLOBAL FUNCTIONS-----------------------------------------*/
extern void forcecheck_nn(Slice *);
extern void setup_BMD(Slice *);
extern void bmdcycle(Slice *psl) ;
/*---------------------------------------------------------------------------*/

/*------------------GLOBAL STRUCTURES-----------------------------------------*/
typedef struct langevin_type {
  // its actually browning dynamics but ok .;) 

    int           step,
                  ninter, //number of langevin steps per propagate_bd
                  bond_breakage_analysis, // 
                  measure_bond_configurations; 

    double        print_time ,
                  total_time,
                  timestep,
                  dtD,
                  friction,
                  diffusionT,
                  diffusionR,
                  mobilityT,
                  sqrtmobilityT,
                  mobilityR,
                  sqrtmobilityR;   // cumulatice time, 
} Langevin;


/*---------------------------------------------------------------------------*/

/*------------------GLOBAL VARIABLES-----------------------------------------*/
Langevin langevin;
/*---------------------------------------------------------------------------*/

#endif
