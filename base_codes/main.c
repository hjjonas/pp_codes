#include "path.h"

Slice        *slice, *psl_old, *start_slice, *copyslice;
System       sys;

vector        nulvec;


// The main function of the simulation program
int main(int argc, char **argv) {

    // Declare variables
    // int icycle, jcycle;     // Loop counters
    // double energy;          // Energy variable
    // FILE *emptyfile, *rdffile;  // File pointers
    time_t t0 = time(0);    // Time at start of program
    
    // Call function to set up simulation parameters
    setup_simulation();
    
    // // Plot the potential if specified in parameters
    // plotpotential();
    
    // // Check the random number generator if specified in parameters
    // check_random();
    
    // // Warm up the Monte Carlo system if specified in parameters
    // if (sys.mc_warmup > 0) {
    //     mc_warmup(&slice[0]);
    // }

    // // Outer loop for the number of cycles specified in parameters
    // for (icycle = 0; icycle < sys.ncycle1; icycle++) {
        
    //     // Print block number
    //     printf("\nBLOCK %d\n", icycle);
        
    //     // Inner loop for the number of cycles specified in parameters
    //     for (jcycle = 0; jcycle < sys.ncycle2; jcycle++) {
            
    //         // If reading path from file is specified in parameters
    //         if (sys.read_path > 0) {
                
    //             // Read coordinates from file
    //             read_coordinates_from_file(&slice[0], sys.read_path*(icycle*sys.ncycle2+jcycle));

    //         // If running Brownian dynamics simulation is specified in parameters
    //         } else if (sys.sim_type == BMD_ALGORITHM) {
                
    //             // Run Brownian dynamics cycle
    //             bmdcycle(&slice[0]);
                
    //             // If printing trajectory is specified in parameters
    //             if (analysis.print_trajectory){
    //                 printing_trajectory(&slice[0]); // try to do it every  1 second
    //             }
            
    //         // If running Monte Carlo simulation is specified in parameters
    //         } else if (sys.sim_type == MC_ALGORITHM) {
                
    //             // Run Monte Carlo cycle
    //             mccycle(&slice[0]);
            
    //         // Otherwise, print error message
    //         } else {
    //             error("choose only MC = sys.sim_type==2 or bmd = sys.sim_type==0"); 
    //         }

    //     }

    //     // Terminate block; printing / measuring 
    //     terminate_block(&slice[0]);

    //     // Print time for this block
    //     time_t block2_1 = time(0);
    //     time_t block2_0 ;
    //     double datetime_diff2 = difftime(block2_1, block2_0);
    //     printf("This block took %lf [sec]\n", datetime_diff2);
    //      block2_0 =block2_1;

    // }

    // // Call function to print final statistics
    // finalstat(&slice[0]);

    // Print total time for program
    time_t t1 = time(0);
    double datetime_diff = difftime(t1, t0);/* in seconds*/
    printf("Total time [min] %lf\n", datetime_diff/60.); 

    // End program
    // free_all_memory();
    return 0;
}

// // this is a very specific function .. 
// void free_all_memory(void){

//     printf("free_all_memory\n");
//     free(&slice[0]);
//     free(&start_slice[0]);
//     freeBondTimeInfoArray(bondtimeinfoarray);

//     if (sys.sim_type==MC_ALGORITHM){
//         free(&psl_old[0]);
//     }

// }

// void bmdcycle(Slice *psl) {
//     /*sys.sim_type==0 Brownina Dynamics*/
//     int istep, update=0, ipart;
//     double dr2;

//     // Perform overdamped langevin dynamics for langevin.step steps
//     for(int istep=0;istep<langevin.step;istep++){
//         propagate_bd(psl);        
//     }

//     // If bond tracking is enabled, perform bond tracking for the current slice
//     if(analysis.bond_tracking==1){
//         bond_tracking(psl);
//         // printBondTimeInfoArray(bondtimeinfoarray )
//     }

//     // Exit the function
// }

// void mc_warmup(Slice *psl){
//  // Use mc_warmup to decorrelate from the input structure; perform MC moves without bond breakages
//     printf("*********  mc_warmup is on:  equilibrating with MC ********\n");

//     // Save old settings to restore later
//     double old_setting=analysis.bond_breakage;
//     double old_beta=psl->beta;
//     int old_clustermc=sys.cluster_MC;
//     int nn=sys.nearest_neighbor;

//     // Set new settings for MC warmup
//     sys.nearest_neighbor=0; // For MC, set to neaghborlist to zero
//     analysis.bond_breakage=0;
//     psl->beta=1.; // Put temperature lower if you want to make stiffer bonds
//     sys.cluster_MC=0;
//     printf("There are %d bonds in the system, keep them fixed as analysis.bond_breakage= %lf.\n",psl->nbonds,analysis.bond_breakage);

//     // Perform MC warmup cycles
//     for(int jcycle=0; jcycle<sys.mc_warmup; jcycle++) {
//         for(int icycle=0; icycle<sys.ncycle2; icycle++) {
//             mccycle(psl); // Perform MC cycle
//         }             

//         // Print status and optimize MC
//         printstatusmc();
//         optimizemc();
//         printenergy_warmup(psl);

//         // Check for bond breakage
//         if (psl->nbonds!=start_slice->nbonds){
//             dprint(psl->nbonds);
//             dprint(start_slice->nbonds);
//             print_slice_information(psl);
//             error("stop. breakage not allowed?");
//         }
//     }

//     // Restore old settings
//     analysis.bond_breakage=old_setting;
//     psl->beta=old_beta; 
//     sys.cluster_MC=old_clustermc;
//     sys.nearest_neighbor=nn;

//     // Update nearest neighbor list if necessary
//     if (sys.nearest_neighbor)  update_nnlist(psl);

//     // Exit the function
//     return;
// }



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
        if (init_conf.restart==0){
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












