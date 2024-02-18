// analysis.h

#ifndef ANALYSIS_H
#define ANALYSIS_H

/*---------------------------------------------------------------------------*/
/*------------------GLOBAL VARIABLES------------------------------------*/


/*---------------------------------------------------------------------------*/
/*------------------GLOBAL FUNCTIONS------------------------------------*/
// extern void  add_bond_information(Slice *,int ,int ,double ,double ,double );
// extern int bond_check(Slice *, int , int );
// extern double particles_distance(Slice *, int , int );
// extern vector particles_vector(Slice *, int , int );
// extern void clustersize_freq_update(Slice *);
// extern void check_maxbonds(Slice *);
// extern int cluster_analysis(Slice *);
// extern void clustersize_identification(Slice *psl);
// extern double running_variance2(double, double, double, double, long);
// extern double running_mean(double, double, long);
// extern void reset_running_statistics(Statistics *);
// extern void running_statistics(Statistics *, double );
// extern int particle_in_wall(Slice *psl, int ipart);
// extern void linked_structure_icluster(Slice *, int , int );
// extern void particles_distance_length_vector(Slice *, int , int  ,double *, double *,vector *);

/*---------------------------------------------------------------------------*/
/*------------------LOCAL FUNCTIONS------------------------------------*/


/*---------------------------------------------------------------------------*/
/*------------------STRUCTURE DEFINITIONS------------------------------------*/

typedef struct statistics_type {
  /* struct for statistics: running mean and variance
  mean = u_n = u_{n-1} + (x_n-u_{n-1})/n
  variance2 = sigma2_n = (sigma2_{n-1} + (x_n-u_n)(x_n-n_{n-1}))/(n-1)*/

    double        mean,
                  variance2;
    long           n;

} Statistics;

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

    Statistics    bin[NPART];

} StatsLength;

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

typedef struct analysis_type {
// list here al the types of analysis the code should do, most scripts are in analysis.c
// and also whether you want to print certain information (print.c)
  int           adjacency,
                s_distribution,
                rdfanalysis,
                bond_op,
                xy_print,
                print_trajectory,
                bond_tracking;
                  
  double        bond_breakage,
                bond_avg;

  int           *s_distribution_ptr,
                *rdfanalysis_ptr,
                *bond_op_ptr;

} Analysis;

/*---------------------------------------------------------------------------*/
/*------------------LOCAL STRUCTURES------------------------------------*/

Cluster cluster;
Analysis analysis;


#endif


