#include "path.h"

Slice        *slice, *cp_slice, *start_slice, *copyslice;
System       sys;
vector       nulvec;



// The main function of the simulation program
int main(int argc, char **argv) {

    // Declare variables
    // int icycle, jcycle;     // Loop counters
    // double energy;          // Energy variable
    // FILE *emptyfile, *rdffile;  // File pointers
    time_t t0 = time(0);    // Time at start of program
    time_t block2_0=t0 ;
    // Call function to set up simulation parameters
    setup_simulation();
    
    // Plot the potential if specified in parameters
    plotpotential();
    
    // Check the random number generator if specified in parameters
    check_random();
    
    // a "warm_up" is performed if you want to attempt to decorrelate the architecture from in input conf. 
    // the bond are not allowed to break or form, but positions/orientations of particles are
    if (init_conf.mc_warmup > 0) {
        mc_warmup(&slice[0]);
    }

    // sys.ncycle1 there are measurements performed (outer loop)
    for (int icycle = 0; icycle < sys.ncycle1; icycle++) {
        
        // Print block number
        printf("\nBLOCK %d\n", icycle);
        
        // Inner loop : decorrelate / propagate the system
        for ( int jcycle = 0; jcycle < sys.ncycle2; jcycle++) {
            switch (sys.sim_type) {
                case READ_TRAJECTORY:
                    read_coordinates_from_file(&slice[0], init_conf.read_path*(icycle*sys.ncycle2+jcycle));
                    break;
                case BMD_ALGORITHM:
                    bmdcycle(&slice[0]);
                    break;
                case MC_ALGORITHM:
                    mccycle(&slice[0]);
                    break;
            }
            // perform version specific analysis here;
            // innerloop_analysis(pls);

        }

        // Terminate block; printing / measuring 
        terminate_block(&slice[0]);

        // Print time for this block
        time_t block2_1 = time(0);
        
        double datetime_diff2 = difftime(block2_1, block2_0);
        printf("This block took %lf [sec]\n", datetime_diff2);
        block2_0 =block2_1;

    }

    // // Call function to print final statistics
    finalstat(&slice[0]);

    // Print total time for program
    time_t t1 = time(0);
    double datetime_diff = difftime(t1, t0);/* in seconds*/
    printf("\n\nTotal time [min] %lf\n", datetime_diff/60.); 

    // End program
    // free_all_memory();
    return 0;
}

void read_coordinates_from_file(Slice *current, int block_number){
    // the file vi trajectroy.xyz  is read
    // each block, separated by a white line (setWHITELINE=1), contains the folloing info in columns:
    // 1: time, 2: particle number, 3-5: position, 6-9: quaternion
    // which are of course NPART lines per block

    static FILE *fp;
    static int n_shots=0;
    static double time_interval=1.;

    int bufferLength=200; // maximum number of characters on line
    char str[bufferLength];
    int ipart,read_lines=0,WHITELINE=0;
    static int linecount=0, init=1;

    double time, x,y,z,time0;
    vector r;
    quaternion q;

    char filename[600];
	snprintf( filename, sizeof( filename ), "%s/trajectory.xyz", init_conf.directorypath);


    // printf("read_coordinates_from_file .. block_number=%d \n",block_number);
    if (init){
    	// open only once in the beginning
	    printf("read_coordinates_from_file %s\n",filename);
    	fp = fopen(filename,"r");
    	init=0;	
    }

    if(fp == NULL) {
      error("Error opening file: trajectory.xyz");
      return;
    }

    // walk over the file line by line
    while(fgets(str, bufferLength, fp)!=NULL ) {
        
        if ((int)(linecount/(current->nparts+WHITELINE)) == block_number){
            // read the block and save it to the slice
            // printf(" reading snapshot number %d \n",block_number);

            // printf("the string: %s",str);
	        // dprint(linecount);
	        // dprint((int)(linecount/(sys.npart )));
	        // dprint(block_number);

            sscanf(str,"%lf %d %lf %lf %lf %lf %lf %lf %lf\n",&time,&ipart,&r.x,&r.y,&r.z,&q.q0,&q.q1,&q.q2,&q.q3);
            pbc(r,sys.boxl);
            current->pts[ipart].r=r;
            current->pts[ipart].q=normalize_quat(q);
            // printf("read: %lf %d %lf %lf %lf %lf %lf %lf %lf\n",time,ipart,x,y,z,q.q0,q.q1,q.q2,q.q3);
            // printf("read: %lf %lf %lf %lf\n",q0,q1,q2,q3);

            //update the patch vectors;
            update_patch_vector_ipart(current,ipart); 
            
            if (read_lines==0){
	            // WE ASSUME THAT  the time_interval  between two frames is 1 second!!
                current->c_time=block_number*time_interval;
            }

            //check for errors
            if (read_lines!=ipart){
                printf("the %dth particle that is read is not ipart=  %d at time=%.1lf \n",read_lines,ipart,current->c_time);
                error("ERROR: there is a problem in read_coordinates_from_file.");     
            }
            
            read_lines++;           
        }
        linecount++;
        if (read_lines==sys.npart){
            // printf("done reading block_number %d \n\n",block_number);
            // reset_center(current); // do this if you want to put com at zero
            // fclose(fp); // don't close the file just yet (its much faster if you only close it at the end)
            n_shots+=1;
            return;
        }
    }

    printf("    you have reached the end of the trajectory file\n");
    if (fgets(str, bufferLength, fp)==NULL){
    	fclose(fp); // close the file now
    	// free_all_memory();
    	error(" end of file of trajectory.");
    }
    
    return;
}

void mc_warmup(Slice *psl){
 // Use mc_warmup to decorrelate from the input structure; perform MC moves without bond breakages
    printf("*********  mc_warmup is on:  equilibrating with MC ********\n");

    // Save old settings to restore later
    int old_setting=analysis.bond_breakage; // why was this a double
    double old_beta=psl->beta;
    int old_clustermc=sys.cluster_MC;
    int nn=sys.nearest_neighbor;
    int xyprint= analysis.xy_print;

    // Set new settings for MC warmup
    sys.nearest_neighbor=0; // For MC, set to neaghborlist to zero
    analysis.bond_breakage=0;
    psl->beta=1.; // Put temperature lower if you want to make stiffer bonds
    sys.cluster_MC=0;
    analysis.xy_print=0;

    printf("There are %d bonds in the system, keep them fixed as analysis.bond_breakage= %d.\n",psl->nbonds,analysis.bond_breakage);

    // Perform MC warmup cycles
    for(int jcycle=0; jcycle<init_conf.mc_warmup; jcycle++) {
        for(int icycle=0; icycle<sys.ncycle2; icycle++) {
            mccycle(psl); // Perform MC cycle
        }             

        // Print status and optimize MC
        printstatusmc();
        optimizemc();
        printenergy(psl, "_warmup");

        // Check for bond breakage
        if (psl->nbonds!=start_slice->nbonds){
            dprint(psl->nbonds);
            dprint(start_slice->nbonds);
            print_slice_information(psl);
            error("stop. breakage not allowed?");
        }
    }

    // Restore old settings
    analysis.bond_breakage=old_setting;
    psl->beta=old_beta; 
    sys.cluster_MC=old_clustermc;
    sys.nearest_neighbor=nn;
    analysis.xy_print=xyprint;
    // Update nearest neighbor list if necessary
    if (sys.nearest_neighbor)  update_nnlist(psl);

    // Exit the function
    return;
}


void emptyfiles(void){
    // deletes files that contain the string ".dat" or, if there is no restart also those files with ".out"

    rmfiles_filename(".", ".dat", 1);
    if (init_conf.restart==0) rmfiles_filename(".", ".out", 1);

    return ;
}

void rmfiles_filename(char directory[MAX_FILENAME_LENGTH], char str_file[MAX_FILENAME_LENGTH], int print){
    struct dirent *de;
    int status;

    if (strlen(directory)==0){  error("rmfile_filename, specify the directory. Current directory is \".\" ");}
    DIR *dr = opendir(directory);

    if (dr == NULL){  // opendir returns NULL if couldn't open directory 
        printf("Could not open current directory. No files are deleted. continue" ); 
        return ; 
    } 
  
    // Refer http://pubs.opengroup.org/onlinepubs/7990989775/xsh/readdir.html 
    // for readdir() 
    while ((de = readdir(dr)) != NULL){
        // printf("%s\n", de->d_name); /* prints all the ".dat" files in the current directory*/
        if (strstr(de->d_name, str_file) != NULL){
            status = remove(de->d_name);

            if (print==1){
                if (status == 0){   printf("deleted successfully:     %s\n", de->d_name);}
                else{               printf("Unable to delete:         %s\n", de->d_name); }
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

    //free memory 
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












