// energy.h

#ifndef ENERGY_H
#define ENERGY_H

/*---------------------------------------------------------------------------*/
/*------------------GLOBAL VARIABLES-----------------------------------------*/

/*---------------------------------------------------------------------------*/
/*------------------STRUCTURE DEFINITIONS------------------------------------*/


/*---------------------------------------------------------------------------*/
/*------------------GLOBAL FUNCTIONS-----------------------------------------*/

extern double total_energy( Slice * );
extern double particle_energy( Slice * , int,int);
extern void energycheck_nn(Slice *);
extern double potential_attractive_bond_energy(Slice *, int , int, double , vector );
extern double gravitational_energy_ipart(Slice *, int );
extern double total_energy_neighborlist(Slice *);
extern double harmonic_oscillator(double , double, double  );
extern void calc_angles( vector ,vector ,vector , double *,double *,double *, int);
extern vector find_directing_patchvector(Slice *, int ,vector );
extern void orientation_parameters(Slice *, int, int , vector, double *, double *, double *);
extern int find_directing_isite(Slice *, int ,vector ); 

/*---------------------------------------------------------------------------*/
/*------------------LOCAL FUNCTIONS------------------------------------------*/
double E_active_alignment(Slice *, int );
void particle_pair_energy(Slice *, int , int , double *, double *, double *);
double total_Epair_neighborlist(Slice *);
double wall_interaction_ipart(Slice *, int );
double total_Epair_woneighborlist(Slice *);

/*---------------------------------------------------------------------------*/
/*------------------GLOBAL STRUCTURES-----------------------------------------*/




#endif
