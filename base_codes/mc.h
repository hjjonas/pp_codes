// mc.h

#ifndef MC_H
#define MC_H

/*------------------STRUCTURE DEFINITIONS------------------------------------*/
typedef struct mcmove_type {
  //these are the different mc move types 
      // each MC type . i.e. single particle or cluster moves rotation or translation , has this structure
      //there are X number of accepter and trial moves, and the ratio = acc/tries 
    long unsigned int           acc,
                                tries;

    double        ratio;

} Mcmove;

typedef struct MC_type {
    // in this structure, you create a MC_type. It tracks the current and final(fin) rot and trans
    // drmax and dqmax are used to optimize the mc moves
    //in name you type e.g. "single particle move" its used when printing output
    char            name[100];

    Mcmove          fintrans,
                    trans,
                    finrot,
                    rot;

    double          drmax, //sigma
                    dqmax; //degrees 

} MC;
/*---------------------------------------------------------------------------*/

/*------------------GLOBAL VARIABLES-----------------------------------------*/
/*---------------------------------------------------------------------------*/

/*------------------GLOBAL FUNCTIONS-----------------------------------------*/
extern void mccycle(Slice *);
extern void setup_MC(void);
extern void optimizemc( );
extern void printstatusmc();
extern void final_printstatusmc_sub(MC );
extern void setup_mc_move(MC *);
/*---------------------------------------------------------------------------*/

/*------------------GLOBAL STRUCTURES-----------------------------------------*/
MC  mc_single, mc_cluster, mc_mono, mc_tailflip;
/*---------------------------------------------------------------------------*/




#endif
