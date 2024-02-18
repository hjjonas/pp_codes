// switchfunction.h

#ifndef SWITCH_H
#define SWITCH_H

/*------------------GLOBAL VARIABLES-----------------------------------------*/
/*---------------------------------------------------------------------------*/

/*------------------STRUCTURE DEFINITIONS------------------------------------*/
typedef struct saccentint_type {
    /*parameters for the integrated Saccent*/
    // swithci_unit used to be switch_variable, watch out with path.inp!
  int             switch_unit;  //0 if the paramters in angles and 1 if in cosangle (only type 1)

  double          switch_expc,          // swich function coefficients (only type 1)
                  switch_expd,
                  switch_expe,
                  switch_expf,
                  switch_expg,
                  switch_exph,
                  switch_expi;

} SaccentInt;


typedef struct site_type {
    //site/swithcfunciton function parameters 

    int           s_accent; // def of S' -> 0 (old switch; broad) 1 (integrated switch, dT dep)                  

    SaccentInt     s_int; // parameters for the integrated patch switch 

     // patch properties ; I think these are used for Sold
    double        S_fixed, // binary patch
                  delta_degree,     /*patch width cutoff in degrees (i.e. where S' is 0)*/
                  cosdelta,         /*patch width cutoff in cosdelta*/
                  oneover_cosdelta,  /*used if s_accent==0 */
                  oneover_cosdeltarev;       /*used if s_accent==0 */  
} Site;
/*---------------------------------------------------------------------------*/

/*------------------GLOBAL FUNCTIONS-----------------------------------------*/
extern double S_value(double , double ,double,int, int);
extern double Saccent(double , int  );
extern double dSaccent_dcosangle(double, int );
extern double dS_dcostheta(Slice *, double , int , double , int  );
/*---------------------------------------------------------------------------*/

/*------------------GLOBAL STRUCTURES-----------------------------------------*/
Site site[NSITES];
/*---------------------------------------------------------------------------*/



#endif
