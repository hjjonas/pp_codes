#include <math.h>
#include <stdint.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

// #include "glut.h"

/*------------------INLINE SUBSTITUTIONS-------------------------------------*/

//Binary definitions
#define TRUE 1
#define FALSE 0
#define CANCELLED -99

#define SIM_MC 2
#define SIM_BMD 0
#define SIM_READ 1


//not so necessary defitions
#define P_INDEX(ptr)  (ptr - pts)
#define F_COMP(a,b) (fabs( a-b ) < 1e-10)
#define MIN(a,b,tmp)    (  (tmp=a) < b  ? tmp : b )
#define MAX(a,b,tmp)    (  (tmp=a) > b  ? tmp : b )
#define MESSAGE(a) printf("Message:"#a"\n")
#define SIGN(a)   ( 2*(a>0) -1) 

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
#define vector_inp(a, b) (a.x * b.x + a.y * b.y + a.z * b.z )
#define vector_cross(a,b,h)   h.x = a.y * b.z - a.z * b.y; h.y = a.z * b.x - a.x * b.z; h.z = a.x * b.y - a.y * b.x;
#define vector_times(a,b,h)   h.x = a.x * b.x; h.y = a.y * b.y; h.z = a.z * b.z;
#define vector_divide(a,b,h)   h.x = a.x / b.x; h.y = a.y / b.y; h.z = a.z / b.z;
#define vector_plustimes(a,b,h)   h.x += a.x * b.x; h.y += a.y * b.y; h.z += a.z * b.z;
#define vector_mintimes(a,b,h)   h.x -= a.x * b.x; h.y -= a.y * b.y; h.z -= a.z * b.z;
#define scalar_times(a,b,h)  h.x = a.x *b; h.y =a.y * b; h.z = a.z * b;  
#define scalar_divide(a,b,h)  h.x = a.x /b; h.y =a.y / b; h.z = a.z / b;  
#define scalar_plustimes(a,b,h)  h.x += a.x *b; h.y +=a.y * b; h.z += a.z * b;  
#define scalar_mintimes(a,b,h)  h.x -= a.x *b; h.y -=a.y * b; h.z -= a.z * b;  
#define vector_add(a,b,h)  {h.x = a.x +b.x; h.y =a.y + b.y; h.z = a.z + b.z;  }
#define vector_minus(a,b,h) { h.x = a.x - b.x; h.y =a.y - b.y; h.z = a.z - b.z;  }
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
    /* never use b and h the same*/
#define quat_add(a,b,h) { h.q0=a.q0+b.q0; h.q1=a.q1+b.q1; h.q2=a.q2+b.q2; h.q3=a.q3+b.q3; }
#define quat_minus(a,b,h) { h.q0=a.q0-b.q0; h.q1=a.q1-b.q1; h.q2=a.q2-b.q2; h.q3=a.q3-b.q3;}
#define quat_inp(a,b) (a.q0*b.q0 + a.q1*b.q1 + a.q2*b.q2 + a.q3*b.q3)
#define sctimes_quat(a,b,h) { h.q0=a.q0*b; h.q1=a.q1*b; h.q2=a.q2*b; h.q3=a.q3*b;}
#define scdivide_quat(a,b,h) { h.q0=a.q0/b; h.q1=a.q1/b; h.q2=a.q2/b; h.q3=a.q3/b; }
#define quat_times(a,b,h) {\
    h.q0 = a.q0*b.q0 - a.q1*b.q1 - a.q2*b.q2 - a.q3*b.q3;\
    h.q1 = a.q0*b.q1 + a.q1*b.q0 + a.q2*b.q3 - a.q3*b.q2;\
    h.q2 = a.q0*b.q2 - a.q1*b.q3 + a.q2*b.q0 + a.q3*b.q1;\
    h.q3 = a.q0*b.q3 + a.q1*b.q2 - a.q2*b.q1 + a.q3*b.q0;} 
#define quat_inverse(a,h) {h.q0 = a.q0; h.q1 = -a.q1; h.q2 = -a.q2; h.q3 = -a.q3;}
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


#define NPART 1000 
#define NSITES 3 //maximum number of sites on 1 particle
#define PTYPES 3  //particle types
#define MAXSLICES 100  // nowhere used? 
#define MAXFRAMES 300 // in bonds breakage code this is 300.000 ; the max number of snapshots to track bond breakage
#define MAXSTATES 3
#define MAXSETS 0
#define NACC 4
#define NSTAT 2
#define MAXBIN 100
#define RDFBINS 100
#define NINTERVALS 250
#define ANGLEBINS 100
#define MAXN0 200
#define MAXK 1000
#define MAXBONDS NSITES
#define MAXTOTBONDS NSITES*NPART



// definitions for neighborlist
#define MAXPARTS            NPART
#define MAX_IMAGES          (MAXPARTS *8) /* is MAXPARTS NPART?*/
#define MAX_NEIGHBORS       (MAXPARTS *300)
#define MAX_IM_NEIGHBORS    300
#define STACKDEPTH          30
#define NSHIFTS             27

#define GRIDOFFSET          3
#define LJ_RZERO            1.5
#define MAXBOXSIZE          20
#define MINCELLSIZE         (LJ_RZERO/GRIDOFFSET)
// #define NCELLX              ((int) (MAXBOXSIZE / MINCELLSIZE)) gives errors??
#define NCELLX              40
#define NCELLY              (NCELLX)
#define NCELLZ              (NCELLX)


/*---------------------------------------------------------------------------*/
/*------------------STRUCTURE DEFINITIONS------------------------------------*/

typedef struct IntArraytype {
  // a list type of intergers
  int *ints;
  size_t used;
  size_t size;

} IntArray ;

typedef struct IntsArraytype {
  // a list type of intergers in 2D
  int *ipart,
      *isite;
  
  size_t used; 
  size_t size; 


} IntsArray ;

typedef struct vector_type {
    double        x,y,z;
} vector;


typedef struct quaternion_type {

    double        q0,q1,q2,q3;

} quaternion;


typedef struct tensor_type {

    vector        x,y,z;

} tensor;

typedef struct statistics_type {
  /* struct for statistics: running mean and variance
  mean = u_n = u_{n-1} + (x_n-u_{n-1})/n
  variance2 = sigma2_n = (sigma2_{n-1} + (x_n-u_n)(x_n-n_{n-1}))/(n-1)*/

    double        mean,
                  variance2;
    long           n;

} Statistics;

typedef struct quattensor_type {

    vector        w,x,y,z;

} quattensor;

typedef struct bondproperty_type{
  /* used for showing nematic order parameter in the graphics for Simons paper
  the bonds are listed with a maximum length of MAXTOTBONDS = NPART*MAXBONDS
  the list contains:
   a location (x,y,z)
   an orientation (q0,q1,q2,q3) in quaterion
   and the nematic order parameter S
  */

  vector       r;

  quaternion   q;

  double       nematic_OP;

} BondProperty;

typedef struct AssDissTimestype {
  /* each patch of each particle has a AssDissTimes struct 
  */

  double  timestamp;    // this is the timestamp of the last association or dissociation 
  
  int     type;    //type is 1 if there was bond formation (association) at timestamp, and 0 for a dissociation path

} AssDissTimes;



typedef struct particle_type {
//info on each particle this pts exists in slice as pts[NPART]
    vector        r,    // position
                  f,    // force BMD
                  t,    // torque BMD
                  dr;   // used in neighborlist 

    quaternion    q;   

    vector        patchvector[NSITES]; //the current ortientation of the patches

    int           cluster_id,
                  ptype, // ptype stands for particle type, with sys.particletype[ptype] you see the characteristics 
                  bonds[MAXBONDS], // saves the particle number of the bond
                  bound_site[MAXBONDS], // saves the site number of ipart that makes the bodn
                  nbonds,
                  a_type2[NSITES],   //associtain dimension 2, dissociation type2 is in bondinfoarray
                  k_bo; //number of bond order properties

    double        bond_energy[MAXBONDS],
                  bond_distance[MAXBONDS];

  
    BondProperty  bond_op[MAXBONDS];

    AssDissTimes  ass_diss[NSITES]; // each patch has a timestamp that indicates its last binding or unbinding 

} Pts;



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


typedef struct activeforce_type {
  /*The active force with strength F0 and direction r [unit vector]*/
  double        F0;

   vector       r;

} ActiveForce;

typedef struct particletype_type {
  /*there are diffent kind of particle type
  e.g. dipatch particle with diameter 3.2 miron = 1 sigma
  or tetrapatch particle with diameter 3.7 micron
  in particle type the number of patches plus diameter (-> gravity) is specified*/

  int           fixed_rot, //fixed rotation 
                fixed_trans, //  fixed   translation
                nparticles,
                nsites,
                activity; //0 no 1 yes (then read in active force)
                

  vector        site[NSITES]; //direction of patch
 
  double        diameter,       /*[sigma]*/
                radius,        /* half of the diameter [sigma]*/
                Rp_sigma,             /*radius of curvature of the patch material [sigma] used for derjaguin scaling */
                delta_rho_kg_m3, /* the density diffence between the lutide-water solution and the colloid [kg/m^3]*/
                fg,           /*gravitational force*/
                zcut,         /*cut off value of z in sigma* = 2mum*/
                b_zc;         /*V_g = LJ(z) or fg*z - b_zc */

  //each particle(type) can only have one active force, it you want different active forces, specify different particles types
  ActiveForce   active_force; 
  
  // each particle(type) has its own association and dissociation rates 
  Statistics    association_times,
                dissociation_times;

} Particletype;

typedef struct bondinfo_type {
  /* the struct that saves the bond infromation, */

  double  s,  //surface-surface distance
          anglei, // angle 1 
          anglej,
          Erep, //repulsive energy
          Ec, // attractive casimir 
          S;   // swithcing function

} BondInfo;





typedef struct statslength_type {
  /* struct for statistics: running mean and variance
  mean = u_n = u_{n-1} + (x_n-u_{n-1})/n
  variance2 = sigma2_n = (sigma2_{n-1} + (x_n-u_n)(x_n-n_{n-1}))/(n-1)*/

    char          filename[100];

    int           length;

    Statistics    *bin;

} StatsLength;

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

typedef struct particles_in_cluster_list_type{
  /* stack,  i.e. a list of particles that are in this cluster, these are not ordered according to the bonds*/
    int           stack[NPART];

} Picl;

typedef struct cluster_type{
  /* the cluster type contain information about all clusters, 
  there are maximally NPART clusters (if all are monomers) 
  in pic (particles in clusterlist) you use the clusterid -> pic[clusterid] 
  to find the stack */

    int           analysis,
                  clustersize[NPART], 
                  update; 

    struct particles_in_cluster_list_type     pic[NPART];

    StatsLength   size_histogram,
                  size_distribution;

} Cluster;


typedef struct slice_type {
// the slice parametrs, in principle you should be able to have mulplie slices
  // with different system settings
    Pts           pts[NPART];

    double        energy,
                  beta,
                  temp,
                  c_time; // current  time 

    int           nclusters,
                  nbonds,
                  nparts;

    Cluster       cluster;
} Slice;



typedef struct langevin_type {
  // its actually browning dynamics but ok .;) 

    int           step,
                  ninter; //number of langevin steps per propagate_bd

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









typedef struct analysis_type {
// list here al the types of analysis the code should do, most scripts are in analysis.c
// and also whether you want to print certain information (print.c)
  int           adjacency,
                s_distribution,
                rdfanalysis,
                bond_op,
                xy_print,
                print_trajectory,
                bond_tracking,
                local_density;
                  
  double        bond_breakage,
                bond_avg;

  int           *s_distribution_ptr,
                *rdfanalysis_ptr,
                *bond_op_ptr;

} Analysis;


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



typedef struct system_type {

    int           particle_setup, 
                  cluster_MC,
                  empty_files,
                  graphics,
                  snapshot,
                  mc_warmup,
                  ncycle1,
                  ncycle2,
                  nearest_neighbor, 
                  npatch,
                  nparticle_types,
                  start_type,
                  sim_type,
                  npart,
                  nchains,
                  restart,
                  read_path,
                  switch_method; //0 is S0, 1 is S90 l; 2 is S lincomb; this holds for all particles;                  

    double        beta,
                  gravity,
                  bond_avg,
                  bond_cutoffE,
                  sigma, /*[mum]*/
                  chaingap;
                  
    vector        boxl;
                  
    Particletype  particletype[PTYPES]; // the type of particles there are

    char          directorypath[500];

} System;

typedef struct cell2_type {
 // for nearest neighbor
  int                n                 ,
                     stack[STACKDEPTH] ;
} Cell;

typedef struct nlist_item_type {
 // for nearest neighbor
  int                ipart             ,
                     image_id          ,
                     *first           ;
} List_item;

typedef struct list_type {
 // for nearest neighbor
  int                cellx             ,
                     celly             ,
                     cellz             ;  

  vector             cellsize          ,
                     inv_cellsize      ,
                     trans[NSHIFTS]    ;

  List_item          nli[MAX_IMAGES]    ;

  int                nlj[MAX_NEIGHBORS] ,
                     nneighbors         ,
                     nimages           ,
                     count             ;

  double             cutoff            ,
                     cutoff2           ,
                     cutoff_buffer2    ,
                     bigpart_cutoff2   ;

  Cell               cell[NCELLX][NCELLY][NCELLZ];

} List;

typedef struct BondTimeInfotype {
  /* BondTimeInfo is a small structure that contains the start time, end time (if officially broken), and the particle numbers
  there is a global (?) integer NActiveBonds that tells you how many "active" bonds there are
  */

  double      start_time,
              end_time;            // timestamp of the timeframe

  int         ipart,
              isite,
              itype2,
              jpart,
              jsite,
              jtype2;

} BondTimeInfo;

typedef struct BondTimeInfoArraytype  {
  // a list type of intergers
  BondTimeInfo *bondinfos;

  size_t used;
  size_t size;

} BondTimeInfoArray ;

typedef struct LocalDensitytype  {
  // the 1D array of the local density  
  // the max density is approx 1 , so local_density \in [0,2 ]     

  StatsLength    local_density;
  
  double         ld_gridsize,
                 box_gridarea; // 

} LocalDensity ;



/*---------------------------------------------------------------------------*/
/*------------------GLOBALLY DEFINED VARIABLES-------------------------------*/


// extern Average rdf_average[RDFBINS], globalCN_average;
// extern Cluster cluster;


extern MC  mc_single, mc_cluster, mc_mono,mc_tailflip;
extern Langevin langevin;
extern List list;

extern Slice *slice,*psl_old,*start_slice,*slices;
extern unsigned int **auto_association_tau0, **auto_dissociation_tau0;
extern StatsLength *auto_association,*auto_dissociation;



extern System sys;
extern vector nulvec,x_axis,y_axis,z_axis;
extern Site site[NSITES]; // used in switchfunction
extern Potential     pot; 
extern Analysis analysis;
// extern BondTimeInfo bondtimeinfo;
extern BondTimeInfoArray *bondinfoarray;
extern IntsArray  *free_sites_ipart;

// extern BondBreakage bondbreakage[MAXFRAMES];
extern LocalDensity ld_array[10];
extern int bound_sites[NPART][NPART],max_tau;
/* PUT THESE IN main.c TOO!*/


/*---------------------------------------------------------------------------*/
/*------------------GLOBALLY DEFINED FUNCTIONS-------------------------------*/


/*version_specific_functions.c*/
extern void local_density(Slice *);
extern void N_TPP_bonds(Slice *);
extern void print_bondtracking(BondTimeInfoArray *,  int  , int , double );
extern void bond_tracking(Slice *);
extern int read_coordinates_from_file(Slice *, int );
extern void identify_bound_and_freesites(Slice *);
extern void print_autocorrelations_tofile(void);
extern void perform_autocorrelation_measurement( int );
extern void update_autocorrelationfunction(StatsLength *, unsigned int *, int );
extern void setup_type2(Slice *);

/*analysis.c*/
extern void  add_bond_information(Slice *,int ,int ,int,double ,double ,double );
extern int bond_check(Slice *, int , int );
extern double particles_distance(Slice *, int , int );
extern vector particles_vector(Slice *, int , int );
extern void clustersize_freq_update(Slice *);
extern void check_maxbonds(Slice *);
extern int cluster_analysis(Slice *);
extern void clustersize_identification(Slice *psl);
extern double running_variance2(double, double, double, double, long);
extern double running_mean(double, double, long);
extern void reset_running_statistics(Statistics *);
extern void running_statistics(Statistics *, double );
extern int particle_in_wall(Slice *psl, int ipart);
extern void linked_structure_icluster(Slice *, int , int );
extern void particles_distance_length_vector(Slice *, int , int  ,double *, double *,vector *);
extern void bond_tracking(Slice *);


/*constructnn.c*/
extern void setup_nnlist();
extern void put_parts_in_box(Slice *);
extern void create_cell_list(Slice *);
extern void find_neighbor_using_cells(Slice *,int);
extern void find_neighbor(Slice *,int );
extern int check_nnlist(Slice *);
extern void update_nnlist(Slice *);
extern void print_nnlist();

/*draw.c*/
extern void reset_center(Slice *);
extern void Init_Graphics();
extern void readslice(Slice *);
extern int readpath();

/*energy.c*/
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

/*force.c*/
extern void gravitational_force_ipart(Slice *, int );
extern void single_bond_force(Slice *, int , int);
extern void propagate_bd(Slice *);
extern void calculate_total_forces(Slice *);
extern void calculate_forces_neighborlist(Slice *);
extern void save_old_positions(Slice *, Slice *);
extern void gravitational_force_ipart(Slice *, int );
extern void forcecheck_nn(Slice *);
extern void derivative_check(Slice *);

/* init.c*/
extern void setup_simulation( void );
extern void init_model(Slice *);
extern void setup_positions_sites( Slice * );
extern void read_particletypes(Slice *);
extern void read_statistics_file(StatsLength *);
extern void plotpotential(void);
extern void check_random(void);

/* main.c*/
extern void terminate_block(Slice *);
extern void finalstat(Slice *);
extern void bmdcycle();
extern void mccycle();
extern void mc_warmup(Slice *);
extern void ffscycle(int, int);
extern void mainloop_for_graphics();
extern void emptyfiles();
extern void error(char *);
extern double cosangle_to_angle(double );
extern int check_nan(double *);
extern void update_patch_vector_ipart(Slice *, int );
extern void update_patch_vectors(Slice *);
extern void free_all_memory(void); 
extern void perform_autocorrelation_measurement(int );
extern void update_autocorrelationfunction(StatsLength *, unsigned int *,int);
extern void allocate_memory(int);

/*mc.c*/
extern void propagate_mc( Slice * );
extern void cluster_propagate_mc(Slice * );
extern void energy_divergence_check(Slice *, char [50]);
extern void optimizemc( );
extern void printstatusmc();
extern int single_particle_move(Slice *, int );
extern int single_particle_move_neighborlist(Slice *, int );
extern int rotatepart_cluster_Ekparts(Slice *);
extern int translatepart_cluster_Ekparts(Slice *);



/*switchfunction.c*/
extern double S_value(double , double ,double,int, int);
extern double Saccent(double , int  );
extern double old_S(double , int);
extern double new_S(double, int);
extern double S_lin(double , int );
extern double S0(double ,double );
extern double S90(double ,double );
extern double S180(double ,double );
extern double switch_method_2( double , double , double , int,int);
extern double dS_dcostheta(Slice *, double , int , double , int  );
extern double dSaccent_dcosangle(double, int );
extern double dSold_dcostheta(double, int  );
extern double dSnew_dcostheta(double , int );
extern double dS0_dcostheta(double, int  ,double, int  );
extern double dS180_dcostheta(double, int  ,double, int  );
extern double dS90_dcostheta(double, int  ,double, int  );



/*terminate_block.c*/
extern void print_input(Slice *);
extern void printrdf();
extern void conf_output(Slice *);
extern void conf_input(Slice *);
extern void print_pos_sites(Slice *);
extern void printenergy(Slice *);
extern void printenergy_warmup(Slice *);
extern void print_statistics_file(StatsLength *, int);
extern void print_statistics_N1(Statistics , char [100]);
extern void print_chain_order(Slice *, int );
extern void print_adjacency_matrix_coordinates(Slice *);
extern void saveScreenshotToFile( int,int);
extern void print_energy_s_theta1_theta2(Slice *);
extern void print_slice_information(Slice *);
extern void printing_trajectory(Slice *);
extern void append_to_file(char [100], int , double );
extern void print_to_file(char [100],int , double);




/*random.c*/
extern double RandomNumber(void);
extern void InitializeRandomNumberGenerator(double);
extern double RandomGaussianNumber();
extern vector RandomBrownianVector(double);
extern vector RandomVector(double);
extern double RandomVelocity(double);
extern vector RandomUnitVector(void);
extern int Random_urandom(void);
extern vector RandomVector3D(vector ) ;
extern double RandomNumberRange(double , double );
extern int RandomIntegerRange(int  , int );
extern vector RandomVector1(void) ;


// criticalcasimir.c
extern double potential_repulsive_energy_sdist(double );
extern void setup_Simons_potential_parameters( void);
extern double potential_repulsive_energy_sdist(double );
extern double bond_repulsive_force(double );
extern double second_der_Vrep(double );
extern double potential_attractive_energy_sdist(double );
extern double bond_attractive_force(double );
extern double second_der_Vc(double );


/*quaternion.c*/
extern double langrange_multiplier_quat(quaternion, quaternion);
extern quaternion quatVecToVec(vector, vector);
extern tensor getrotmatrix(quaternion);
extern quattensor getquatmatrix(quaternion);
extern quaternion RandomQuaternion();
extern quaternion RandomQuaternionRange(double);
extern quaternion QuaternionXaxis( double );
extern quaternion QuaternionYaxis( double );
extern quaternion QuaternionZaxis( double );
extern void rotate_quaternion_ipart( Slice *, int , quaternion  );
extern quaternion FlipQuaternion(vector  , double ) ;
extern quaternion normalize_quat(quaternion );

/*storage.c*/
extern void memory_allocation(void);
extern void initIntArray(IntArray *, size_t );
extern void insertIntArray(IntArray *, int ); 
extern void freeIntArray(IntArray *);
extern void removeIthElementIntArray(IntArray *, int  );
extern void removeElementXIntArray(IntArray *, int  );
extern int checkElementXIntArray(IntArray *, int );
extern void printIntArray(IntArray * );

extern void initIntsArray(IntsArray *, size_t );
extern void insertIntsArray(IntsArray *, int ,int  ); 
extern void freeIntsArray(IntsArray *);
extern void removeIthElementIntsArray(IntsArray *, int  );
extern void removeElementXIntsArray(IntsArray *, int, int  );
extern int checkElementXIntsArray(IntsArray *, int , int );
extern void printIntsArray(IntsArray * );

extern void initBondTimeInfoArray(BondTimeInfoArray *, size_t );
extern void insertBondTimeInfoArray(BondTimeInfoArray *, BondTimeInfo ); 
extern void freeBondTimeInfoArray(BondTimeInfoArray *);
extern void removeIthElementBondTimeInfoArray(BondTimeInfoArray *, int  );
extern void removeElementXBondTimeInfoArray(BondTimeInfoArray *, BondTimeInfo  );
extern int checkElementXBondTimeInfoArray(BondTimeInfoArray *, int ,int );
extern void memory_allocation(void);
extern void printBondTimeInfoArray(BondTimeInfoArray * );
extern int checkElementXBondTimeInfoArray_atIth(BondTimeInfoArray *, int ,int , int);
extern void insertBondTimeInfoArray_atIth(BondTimeInfoArray *, BondTimeInfo ,int); 
