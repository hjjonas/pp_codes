// criticalcasimir.h

#ifndef CC_H
#define CC_H


/*---------------------------------------------------------------------------*/
/*------------------GLOBALLY DEFINED FUNCTIONS-------------------------------*/
extern double potential_repulsive_energy_sdist(double );
extern void setup_criticalCasimir_potential_parameters( void);
extern double potential_repulsive_energy_sdist(double );
extern double bond_repulsive_force(double );
extern double second_der_Vrep(double );
extern double potential_attractive_energy_sdist(double );
extern double bond_attractive_force(double );
extern double second_der_Vc(double );

/*---------------------------------------------------------------------------*/
/*------------------STRUCTURE DEFINITIONS------------------------------------*/

typedef struct potential_type {
    //critical Casimir parameters

    double        r_wetting,
                  surface_charge,
                  dT,
                  A,
                  xi,
                  Arep,
                  dl;

    //potential properties
    double        s_overlap, // the distance at which the system overlaps. Different for casimir and LJ
                  s_min, // s at minimum energy V=Vrep+Vattr
                  s_forcecut,
                  epsilongravLJ, // for gravity
                  S_fixed,
                  Erep_smin,
                  Ec_smin,
                  E_smin,
                  bond_cutoffE,
                  s_cutoff,
                  rcutoff;

    // gravity potential parameters
    double        sqw_epsilon,
                  wall_int;  

} Potential;


/*---------------------------------------------------------------------------*/
/*------------------LOCAL STRUCTURES------------------------------------*/

Potential pot;



#endif