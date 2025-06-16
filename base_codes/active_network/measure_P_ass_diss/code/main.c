#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "path.h"
#include <stdint.h>
#include <time.h>
#include <unistd.h>
#include <dirent.h> 
// #include <omp.h> // for paralelization


/******* THIS CODE IS TO READ TRAJECTORIES FROM FILE AND ANALYSE THEM ******/
//  for example, measure the local density, autocorrelation of the association and dissocition, and bond lifetime

LocalDensity ld_array[10];
MC  mc_single, mc_cluster, mc_mono,mc_tailflip;
Langevin langevin;
List list;

Slice *slice,*psl_old, *start_slice, *slices;

unsigned int **auto_association_tau0, **auto_dissociation_tau0;
StatsLength *auto_association,*auto_dissociation;



int max_tau;
System sys;
vector nulvec;
Site site[NSITES];
Potential     pot;
Analysis analysis;
BondTimeInfoArray *bondinfoarray;
IntsArray  *free_sites_ipart;

int bound_sites[NPART][NPART]={-1};

// In this main function, we read the data and analyze the data
int main(int argc, char **argv) {

    time_t t0 = time(0);    // Time at start of program
    
    // initialization of simulation 
    setup_simulation();
    
    // Plot the potential to see its shape 
    plotpotential();
    
    if (sys.read_path <= 0) {
        error("sys.read_path must be > 0");  
    }

    // start_t is the first timestamp you want to read.
    //  skip is set by   read_path and is in 
    int check,tau_start,tau,max_tau_old,start_t=0;
    max_tau=10000; // my maximum length of (physicsal) time is 10000*1 seconds 
    max_tau_old=max_tau;

    printf(" first load in all data  from file \n");
    slices = (Slice *)calloc(max_tau,sizeof(Slice));
    int skip=sys.read_path;
    printf(" step 1.  starting at timestamp %d, timestep sizes of %d seconds , and #frames=%d \n",start_t,skip ,max_tau);
    dprint(skip);
    for (tau=0;tau<max_tau;tau++){ 
        check=read_coordinates_from_file(&slices[tau], start_t+skip*tau);
        slices[tau].energy=total_energy(&slices[tau]);

        if (check==-1){ // -1 is returned if the end-of-file has been reached. then save max_tau as tau-1
            max_tau=(int)(tau-1);
            break;
        }
        
        if (analysis.local_density) local_density(&slices[tau]);
    }

    printf("      done reading the data. start_t=%d,max_tau=%d \n",start_t,max_tau);
    allocate_memory(max_tau);

    printf("    step 2. perform_autocorrelation_measurement \n");   
    perform_autocorrelation_measurement( max_tau);
    print_autocorrelations_tofile();
    
    // Print total time for code to run
    time_t t1 = time(0);
    double datetime_diff = difftime(t1, t0);/* in seconds*/
    printf("Total time [min] %lf\n", datetime_diff/60.); 

   
    free_all_memory();
    return 0;
}

void allocate_memory(int max_tau){

    int reaction_types=4;

    // there are 4 types of binding/unbinding  reactions 
    // 1) DP + DP -> DP-DP
    // 2) TPP* + DP -> TPP*-DP   here the active force points TOWARD the bond formed
    // 3) TPP* + DP -> *TPP-DP   here the active force points AWAY   the bond formed
    // 4) all of the above together


    // Allocate memory for auto_association and auto_dissociation arrays
    auto_association_tau0 = (unsigned int **)calloc(reaction_types , sizeof(unsigned int *));
    auto_dissociation_tau0 = (unsigned int **)calloc(reaction_types , sizeof(unsigned int *));

     // Allocate memory for bin arrays inside auto_association and auto_dissociation
    auto_dissociation = (StatsLength *)calloc(reaction_types, sizeof(StatsLength));
    auto_association = (StatsLength *)calloc(reaction_types, sizeof(StatsLength));


    for (int n = 0; n < reaction_types; n++) {
        
        auto_association_tau0[n] = (unsigned int *)calloc(max_tau, sizeof(unsigned int));
        auto_dissociation_tau0[n] = (unsigned int *)calloc(max_tau, sizeof(unsigned int));

        auto_dissociation[n].bin = (Statistics *)calloc(max_tau, sizeof(Statistics));
        auto_association[n].bin = (Statistics *)calloc(max_tau, sizeof(Statistics));

        auto_dissociation[n].length = max_tau;
        auto_association[n].length = max_tau;
        
    }
    // printf("done\n");

}
void free_all_memory(void){

    free(&slice[0]);
    free(&start_slice[0]);
    freeBondTimeInfoArray(bondinfoarray);
    int reaction_types=4;
    for (int n = 0; n <reaction_types;n++) {
        free(auto_association[n].bin);
        free(auto_dissociation[n].bin);

        free(auto_association_tau0[n]);
        free(auto_dissociation_tau0[n]);    
    }

    // Free the arrays of pointers
    free(auto_association_tau0);
    free(auto_dissociation_tau0);
    free(auto_association);
    free(auto_dissociation);

    // freeIntArray(free_sites_isite) ;
    freeIntsArray(free_sites_ipart) ;
    free(slices);

    

    if (sys.sim_type==SIM_MC){
        free(&psl_old[0]);
    }

}



void emptyfiles(){
    struct dirent *de;
    int status;

    DIR *dr = opendir(".");
    if (dr == NULL)  // opendir returns NULL if couldn't open directory 
    { 
        printf("Could not open current directory. No files are deleted. continue" ); 
        return ; 
    } 
  
    // Refer http://pubs.opengroup.org/onlinepubs/7990989775/xsh/readdir.html 
    // for readdir() 
    while ((de = readdir(dr)) != NULL){
        // printf("%s\n", de->d_name); /* prints all the ".dat" files in the current directory*/
        if (strstr(de->d_name, ".dat") != NULL){
            status = remove(de->d_name);

            if (status == 0){
                printf("deleted successfully:     %s\n", de->d_name);
            }
            else{
                printf("Unable to delete the file\n");
            }
        }
        // if you want to don't restart from a time and configuration, also remove ".out" files
        if (sys.restart==0){
            if (strstr(de->d_name, ".out") != NULL){
                status = remove(de->d_name);
    
                if (status == 0){
                    printf("deleted successfully:     %s\n", de->d_name);
                }
                else{
                    printf("Unable to delete the file\n");
                }
            }
        }
    }
  
    closedir(dr);     
    return ;
}

void error(char *msg){
/*This function is a simple error-handling function that takes an error message as an argument 
    and prints it along with the "error:" prefix. 
    It then exits the program with an error code of 0. */
  printf("error: %s\n",msg);

  // free(&slice[0]);
  // free(&start_slice[0]);

  free_all_memory();

  exit(0);
}


double cosangle_to_angle(double cosangle){
    // acos(cosangle) gives a nan if cosangle>1
    
    double angle, rad_to_degree=180./PI;
    angle=(cosangle<1.)? acos(cosangle)*rad_to_degree:0;
    if (cosangle<=-1){
        angle=180;
    }
    // return  the angle in degrees
    return angle;
}



int check_nan(double *number){
/*This function checks if a given pointer points to a NaN (not-a-number) value. 
    If so, it prints the value and the name of the variable,
    then calls the error() function to terminate the program with an error code of 0. */
    if (isnan(*number)!=0){
        gprint(*number);
        printf("%s:\n",getName(*number));

        error("this number is nan ");
    }
    return 0;
}


void update_patch_vector_ipart(Slice *psl, int ipart){
    // reads the quaternion of the particle and transforms it into a patch vector for each site
    Pts *psi;
    int isite,ptype;
    tensor rotmati;
    
    psi = &psl->pts[ipart]; // particle ipart
    rotmati = getrotmatrix(psi->q);  // return the rotation matrix based on the quaterion
    for( isite=0; isite<sys.particletype[psi->ptype].nsites; isite++) {
        ptype=psi->ptype;
        // in sys.particletype[psi->ptype].site the standard patch vector is written.
        matrix_x_vector(rotmati,sys.particletype[psi->ptype].site[isite],psi->patchvector[isite]); 
    }
    return;
}

void update_patch_vectors(Slice *psl){
    // in the quaternion is updates, you also need to update the patch vectors, 
    // this funcitons updates all patch vectors
    int ipart;
    for(ipart=0;ipart<psl->nparts;ipart++){
        update_patch_vector_ipart(psl, ipart);
    }
    return;
}












