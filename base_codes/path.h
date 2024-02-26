#include <math.h>
#include <stdint.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <dirent.h> 
#include <stdbool.h>
#include <time.h>

/*--------------------------------------------------------------------------------------*/
/*------------------FIXED/HARD CODED MAXIMUM VALUES ------------------------------------*/

#define NPART 1000        // maximum number of particles
#define NSITES 6          // maximum number of sites on 1 particle
#define PTYPES 6          // particle types; eg dipatch and tripatch cC potential 
#define MAXBONDS NSITES   // maximum number of bonds per particle

// to select way to "move the particles" with sys.sim_type
#define BMD_ALGORITHM 1  
#define MC_ALGORITHM 2
#define READ_TRAJECTORY 3


/*---------------------------------------------------------------------------*/
/*------------------STRUCTURE DEFINITIONS------------------------------------*/

typedef struct vector_type { 
    double        x,y,z;
} vector;

typedef struct quaternion_type {
  // see section 2.4.3 RigidBodyRotations&Quaternions of thesis HJJ
    double        q0,q1,q2,q3;
} quaternion;

typedef struct particletype_type {
  /*there are diffent kind of particle type
  e.g. dipatch particle with diameter 3.2 miron = 1 sigma
  or tetrapatch particle with diameter 3.7 micron
  in particle type the number of patches plus diameter (-> gravity) is specified*/

  int           nparticles,
                nsites,
                activity; //0 no 1 yes (then read in active force)
                
  vector        site[NSITES], //direction of patch          
                e_A;          // active force direction, note: each particle(type) can only have one active force, it you want different active forces, specify different particles types

  double        diameter,       /*[sigma]*/
                radius,        /* half of the diameter [sigma]*/
                Rp_sigma,             /*radius of curvature of the patch material [sigma] used for derjaguin scaling */
                delta_rho_kg_m3, /* the density diffence between the lutide-water solution and the colloid [kg/m^3]*/
                fg,           /*gravitational force*/
                zcut,         /*cut off value of z in sigma* = 2mum*/
                b_zc,         /*V_g = LJ(z) or fg*z - b_zc */
                F_A;          // active force magnitude
  
  
  // each particle(type) has its own association and dissociation rates 
  // Statistics    association_times,
  //               dissociation_times;

} Particletype;

typedef struct particle_type {
//info on each particle this pts exists in slice as pts[NPART]
    vector        r,    // position
                  f,    // force BMD
                  t,    // torque BMD
                  dr;   // used in neighborlist, it tracks how far the particlce diffused wrt moment at which neighborlist is created

    quaternion    q;    // the current quaterion of the particle

    vector        patchvector[NSITES];  // the current ortientation of the patches; note you need to apply "INSERT FUNCTION" after new q

    int           cluster_id,           // the id of the cluster, cluster= all particles that are connected by bonds. see in energy.c what the definition of a bond is 
                  ptype,                // ptype stands for particle type, with sys.particletype[ptype] you see the characteristics 
                  bonds[MAXBONDS],      // here the particlenumbers of the bound particles are listed
                  nbonds;               // the number of bonds 
                  // k_bo;                 // number of bond order properties, not used? 

    double        bond_energy[MAXBONDS],    // the bond energy of the ith bond
                  bond_distance[MAXBONDS];  // the r_ss of the ith bond

    // specific version of active_network, add them later to versions of the code
    // BondProperty  bond_op[MAXBONDS];
    // AssDissTimes  ass_diss[NSITES]; // each patch has a timestamp that indicates its last binding or unbinding 

} Pts;

typedef struct system_type {

    int           cluster_MC, // PUT IN MC STRUCT?
                  ncycle1,
                  ncycle2,
                  nearest_neighbor,  // put it in BMD?? not used in MC
                  npatch,           // in particletype?
                  nparticle_types,
                  sim_type,       // Identify here if you use MC or BMD or something else
                  npart,
                  switch_method;  //  0 is S0, 1 is S90 l; 2 is S lincomb; this holds for all particles;                  

    double        beta,           // inverse temperature beta
                  gravity,        // on (1) or off (0)
                  bond_avg,       // REPLACE, this is a measurement, not a system variable
                  bond_cutoffE,   // the energy that defines the energetic cut-off for defining a bond 
                  sigma;          /* [mum]*/

    vector        boxl;                 // the (periodic) box size (Lx,Ly,Lz)
                  
    Particletype  particletype[PTYPES]; // the type of particles there are

} System;

typedef struct slice_type {
// the slice parametrs, in principle you should be able to have mulplie slices
  // with different system settings

    Pts           pts[NPART];

    double        energy,
                  beta,
                  temp,
                  bond_probability,   // this is a measurement on the slice, maybe make separate struct for the measurements
                  c_time; // current  time 

    int           nclusters,
                  nbonds,
                  nparts;

} Slice;

/*---------------------------------------------------------------------------*/
/*------------------THE H-FILES----------------------------------------------*/

#include "criticalcasimir.h"
#include "storage.h"
#include "analysis.h"
#include "constructnn.h"
#include "init.h"
#include "random.h"
#include "quaternion.h"
#include "switchfunction.h"
#include "terminate_block.h"
#include "mc.h"
#include "energy.h"
#include "bmd.h"
#include "version_specific_functions.h"

/*---------------------------------------------------------------------------*/
/*------------------GLOBALLY DEFINED VARIABLES-------------------------------*/

extern Slice        *slice, *psl_old, *start_slice, *copyslice;
extern System       sys;
extern vector       nulvec;

// in criticalcasimir.h
extern Potential     pot; 

// in analysis.h
extern Cluster cluster;
extern Analysis analysis;

// in quaternion.h

// in mc.h
extern MC  mc_single, mc_cluster, mc_mono, mc_tailflip;

// in switchfunction.h
extern Site site[NSITES]; 

// in constructnn.h
extern List list;

// in bmd.h
extern Langevin langevin;

/*---------------------------------------------------------------------------*/
/*------------------GLOBALLY DEFINED FUNCTIONS-------------------------------*/

/* main.c*/
// extern void terminate_block(Slice *);
// extern void finalstat(Slice *);
// extern void bmdcycle();
extern void mc_warmup(Slice *);
// extern void ffscycle(int, int);
// extern void mainloop_for_graphics();
extern void emptyfiles();
extern void error(char *);
extern double cosangle_to_angle(double );
extern int check_nan(double *);
extern void update_patch_vector_ipart(Slice *, int );
extern void update_patch_vectors(Slice *);
extern void read_coordinates_from_file(Slice *, int );
extern void rmfiles_filename(char [MAX_FILENAME_LENGTH], char [MAX_FILENAME_LENGTH], int );
// extern void free_all_memory(void); 



/*------------------INLINE SUBSTITUTIONS-------------------------------------*/

//Binary definitions
#define TRUE 1
#define FALSE 0
#define CANCELLED -99

#define P_INDEX(ptr)    (ptr - pts)
#define F_COMP(a,b)     (fabs( a-b ) < 1e-10)
#define MIN(a,b,tmp)    (  (tmp=a) < b  ? tmp : b )
#define MAX(a,b,tmp)    (  (tmp=a) > b  ? tmp : b )
#define MESSAGE(a)      printf("Message:"#a"\n")
#define SIGN(a)         ( 2*(a>0) -1) 

#define max(a,b) \
    ({ __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a > _b ? _a : _b; })
#define min(a,b) \
    ({ __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a < _b ? _a : _b; })

//Math shorthands + special numbers
#define PI           3.141592653589793
#define HALFPI       1.5707963267948966
#define SQRT2        1.4142136
#define BIGNUM       1e99
#define EPS          1e-10
#define NVT          1
#define cubic(x) ((x)*(x)*(x))
#define square(x) ((x)*(x))

//Print shorthands
#define dprint(expr) printf("  "#expr " = %d\n",expr)
#define gprint(expr) printf("  "#expr " = %20.15g\n",expr)
#define vprint(expr) printf("  "#expr " = ( %16.14g %16.14g %16.14g ) \n",expr.x, expr.y, expr.z)
#define qprint(expr) printf("  "#expr " = ( %16.14g %16.14g %16.14g %16.14g) \n",expr.q0, expr.q1, expr.q2 expr.q3)
#define getName(var)  #var

//Vector definitions
#define vector_inp(a, b)        (a.x * b.x + a.y * b.y + a.z * b.z )
#define vector_cross(a,b,h)     h.x = a.y * b.z - a.z * b.y; h.y = a.z * b.x - a.x * b.z; h.z = a.x * b.y - a.y * b.x;
#define vector_times(a,b,h)     h.x = a.x * b.x; h.y = a.y * b.y; h.z = a.z * b.z;
#define vector_divide(a,b,h)    h.x = a.x / b.x; h.y = a.y / b.y; h.z = a.z / b.z;
#define vector_plustimes(a,b,h) h.x += a.x * b.x; h.y += a.y * b.y; h.z += a.z * b.z;
#define vector_mintimes(a,b,h)  h.x -= a.x * b.x; h.y -= a.y * b.y; h.z -= a.z * b.z;
#define scalar_times(a,b,h)     h.x = a.x *b; h.y =a.y * b; h.z = a.z * b;  
#define scalar_divide(a,b,h)    h.x = a.x /b; h.y =a.y / b; h.z = a.z / b;  
#define scalar_plustimes(a,b,h) h.x += a.x *b; h.y +=a.y * b; h.z += a.z * b;  
#define scalar_mintimes(a,b,h)  h.x -= a.x *b; h.y -=a.y * b; h.z -= a.z * b;  
#define vector_add(a,b,h)       {h.x = a.x +b.x; h.y =a.y + b.y; h.z = a.z + b.z;  }
#define vector_minus(a,b,h)     { h.x = a.x - b.x; h.y =a.y - b.y; h.z = a.z - b.z;  }
#define matrix_x_vector(m,a,h) {\
    h.x = m.x.x*a.x + m.x.y*a.y + m.x.z*a.z;\
    h.y = m.y.x*a.x + m.y.y*a.y + m.y.z*a.z;\
    h.z = m.z.x*a.x + m.z.y*a.y + m.z.z*a.z;}
#define matrixT_x_vector(m,a,h) {\
    h.x = m.x.x*a.x + m.y.x*a.y + m.z.x*a.z;\
    h.y = m.x.y*a.x + m.y.y*a.y + m.z.y*a.z;\
    h.z = m.x.z*a.x + m.y.z*a.y + m.z.z*a.z;}
#define normvec(a,b) {b.x=a.x/sqrt(vector_inp(a,a)); b.y=a.y/sqrt(vector_inp(a,a)); b.z=a.z/sqrt(vector_inp(a,a));}

//Quaternion definitions
    /* never use b and h the same variable*/
#define quat_add(a,b,h)     { h.q0=a.q0+b.q0; h.q1=a.q1+b.q1; h.q2=a.q2+b.q2; h.q3=a.q3+b.q3; }
#define quat_minus(a,b,h)   { h.q0=a.q0-b.q0; h.q1=a.q1-b.q1; h.q2=a.q2-b.q2; h.q3=a.q3-b.q3;}
#define quat_inp(a,b)       (a.q0*b.q0 + a.q1*b.q1 + a.q2*b.q2 + a.q3*b.q3)
#define sctimes_quat(a,b,h) { h.q0=a.q0*b; h.q1=a.q1*b; h.q2=a.q2*b; h.q3=a.q3*b;}
#define scdivide_quat(a,b,h){ h.q0=a.q0/b; h.q1=a.q1/b; h.q2=a.q2/b; h.q3=a.q3/b; }
#define quat_times(a,b,h) {\
    h.q0 = a.q0*b.q0 - a.q1*b.q1 - a.q2*b.q2 - a.q3*b.q3;\
    h.q1 = a.q0*b.q1 + a.q1*b.q0 + a.q2*b.q3 - a.q3*b.q2;\
    h.q2 = a.q0*b.q2 - a.q1*b.q3 + a.q2*b.q0 + a.q3*b.q1;\
    h.q3 = a.q0*b.q3 + a.q1*b.q2 - a.q2*b.q1 + a.q3*b.q0;} 
#define quat_inverse(a,h)   {h.q0 = a.q0; h.q1 = -a.q1; h.q2 = -a.q2; h.q3 = -a.q3;}
#define quatmatrix_x_vec(m,a,h) {\
    h.q0 = m.w.x*a.x + m.w.y*a.y + m.w.z*a.z;\
    h.q1 = m.x.x*a.x + m.x.y*a.y + m.x.z*a.z;\
    h.q2 = m.y.x*a.x + m.y.y*a.y + m.y.z*a.z;\
    h.q3 = m.z.x*a.x + m.z.y*a.y + m.z.z*a.z;}

//Update definitions
#define update_average(a,b) {a.now =b; a.sum += a.now; a.sumsq += a.now*a.now; a.n++;}
#define update_blockaver(a) {a.sum /= a.n; a.sumsq = sqrt(a.sumsq/a.n - a.sum*a.sum);}
#define update_finaver(a,b) {a.sum += b.sum; a.sumsq += b.sum*b.sum;  b.sum = 0; b.sumsq = 0; b.now=0; b.n=0;}
#define update_blockmc(a,b) {a.acc += b.acc; a.tries += b.tries;  b.acc = 0; b.tries = 0;}


#define SWAP(a,b,c) {c=a;a=b;b=c;}

//Periodic Boundary conditions
#define pbc(a,b) {\
    do{\
      if(a.x>(0.5*b.x)) a.x-=b.x;\
      if(a.x<(-0.5*b.x)) a.x+=b.x;\
    }while( fabs(a.x)>0.5*fabs(b.x) );\
    do{\
      if(a.y>(0.5*b.y)) a.y-=b.y;\
      if(a.y<(-0.5*b.y)) a.y+=b.y;\
    }while( fabs(a.y)>0.5*fabs(b.y) );\
    do{\
      if(a.z>(0.5*b.z)) a.z-=b.z;\
      if(a.z<(-0.5*b.z)) a.z+=b.z;\
      }while( fabs(a.z)>0.5*fabs(b.z) );\
    }
