#include "path.h"

// ***************  all function specific for simulating the cluster size distribution or chain configurations ***************

/*------------------LOCAL FUNCTIONS------------------------------------------*/
void print_xy_positions(Slice *); // for the chains 
/*---------------------------------------------------------------------------*/



void version_specific_analysis(Slice *psl){
	// this funciton is located inside terminate_block() at the outerloop in main.c	

    clustersize_freq_update(psl);
	
	if(cluster.analysis==1){
        // the histogram is the same as size_distribution, but not normalized
        print_StatsLength_to_file(&cluster.size_distribution);
    }

    // set to 1 if you want to print the xy postion of the chain to file
    //  make sure that analysis.bond_breakage=0 to keep the chain intact
    if ((analysis.xy_print==1) && (analysis.bond_breakage==0)) print_xy_positions(psl);

}

void innerloop_analysis(Slice *psl){
	// this funciton is located inside in innerloop in main.c
	// so after ncycle2 bmd, mc or trajectory reads, you perform this function. 
}

void clustersize_freq_update(Slice *psl){
    /*this function saves the frequency of occurence of the clusterlengths  */

    int i,l;
    int   freq_clusterlength[NPART]={0};
    
    for(i=0;i<psl->nclusters;i++){
        freq_clusterlength[cluster.clustersizes[i]-1]++; // use -1 here because size one is the smallest size
    } 

    for(l=0;l<psl->nparts;l++){
        /*loops over the lengths which has maximumvalue of psl->nparts*/
        // update_average(chain.length_histogram[i], freq_clusterlength[i]);
        running_statistics(&cluster.size_histogram.bin[l], freq_clusterlength[l]);
        running_statistics(&cluster.size_distribution.bin[l], (double)freq_clusterlength[l]/psl->nclusters);
    } 

    return; 
}


void print_xy_positions(Slice *psl){
    // x   y   fx  fy  particle    time
    int ipart;
    static int linecount_xypos=0;
    static double cycle_mc=0;
    FILE *posfile;
    char filename[100];
    char* a = "xypos";
    char* extension = ".csv";

    vector v1,r_ref,dv,dv_norm,dr;
    vector yaxis=nulvec,xaxis=nulvec,drcheck;
    vector v[psl->nparts],v_new,f_new;
    double costheta,rotangle;
    quaternion q;
    tensor R;
    vector rotvec;

    yaxis.y=1.;
    xaxis.x=1.;


    snprintf( filename, sizeof( filename ), "%s%s", a,  extension );
    
    if( access( filename, F_OK ) == -1 ) {
        posfile = fopen(filename,"w");
        if (posfile == NULL){
            printf("Error with creating \"%s\" file",filename);
        }
        else{
            printf("creating the file %s",filename);
            fprintf(posfile,",x,y,fx,fy,particle,time\n");
        }
        fclose(posfile);
    }
    
    posfile = fopen(filename,"a");
    if (posfile == NULL){
        printf("Error with opening \"%s\" file",filename);
    }
    else{
        printf("\nPrinting the xy positions to %s\n",filename);
        /* rotate the chain before printing, such that the chain is oriented along the x-axis. This is better for the analysis tool of Simon*/
        /*  translate all particles such that v0=(0,0,0) 
            take the angle between  the x-value of particle 0 along x-axis; the y-value should be around 0.
            and the vector between particle 0 and N-1
            use this angle to do the rotation via the quaternion in z-axis QuaternionZaxis(degree) */
       
        // calculate ther rotation matrix R
        
        // shift the reference frame. walk through the chain
        r_ref=psl->pts[0].r;
        v[0]=nulvec;// /choose particle 0 as reference
        // vprint(v[0]);
        for(ipart=1; ipart<psl->nparts; ipart++) {
            //the vector difference between particle i and j
            vector_minus(psl->pts[ipart].r,psl->pts[ipart-1].r,dr);
            pbc(dr,sys.boxl);
            vector_add(v[ipart-1],dr,v[ipart]);
            // vprint(v[ipart]);
        }
        
        // calculate the rotation matrix
        dv.x=v[psl->nparts-1].x;
        dv.y=v[psl->nparts-1].y;
        dv.z=0.;

        scalar_divide(dv,sqrt(vector_inp(dv,dv)),dv_norm); 

        /*build in if angle>180(== costheta<0, the rotation becomes is 360-angle)*/
        costheta=vector_inp(dv_norm,xaxis);
        if (dv_norm.y>0.){
            //rotations are done tegen de klok in
            rotangle=360.-acos(costheta)*180./PI;

        }
        else{
            rotangle=acos(costheta)*180./PI;
        }
        // gprint(rotangle);
        q=QuaternionZaxis(rotangle);
        R = getrotmatrix(q);
        
        //perform the rotarion
        for(ipart=0; ipart<psl->nparts; ipart++) {
            matrix_x_vector(R, v[ipart], v_new);

            //print to file        
            if (sys.sim_type==MC_ALGORITHM){ //mc
                fprintf(posfile,"%d,%.5lf,%.5lf,0.0,0.0,%d,%.5lf\n", linecount_xypos,v_new.x, v_new.y,ipart ,cycle_mc);
                linecount_xypos++;
            }   
            else if(sys.sim_type==BMD_ALGORITHM){ //bmd
                matrix_x_vector(R, psl->pts[ipart].f, f_new);
                // vprint(psl->pts[ipart].f);
                // gprint(vector_inp(psl->pts[ipart].f,psl->pts[ipart].f));
                // vprint(f_new);
                // gprint(vector_inp(f_new,f_new));
                fprintf(posfile,"%d,%.5lf,%.5lf,%.5lf,%.5lf,%d,%.5lf\n", linecount_xypos ,v_new.x, v_new.y,f_new.x,f_new.y,ipart ,psl->c_time);
                linecount_xypos++;
            }         
        }
    } 
    
    cycle_mc+=0.1;
    fclose(posfile);
    return;

}



/*_________________________MEMORY ALLOCATION AND FREEING________________*/

void special_init(Slice *psl){
    // for cluster analysis 
    if(cluster.analysis==1){
        // chain. = (Statistics *)calloc(NPART,sizeof(Statistics));
        printf("\nThe code performs cluster analysis \n");
        printf("      the total energy is %lf\n", total_energy(psl));
        printf("      now cluster analysis\n");
        psl->nclusters = cluster_analysis(psl); 
        clustersize_identification(psl);

        printf("      Initial # bonds is %d \n", psl->nbonds);
     
        strcpy(cluster.size_histogram.filename,"clustersize_histogram.out");
        strcpy(cluster.size_distribution.filename,"clustersize_distribution.out");
 
    }
    
}
void allocate_memory(Slice *psl){
	// allocate versions specific memory for cluster anallysys 

	// dynamically allocate memory for cluster.size_distribution.bin 
	cluster.size_histogram.bin = (Statistics *)calloc(psl->nparts, sizeof(Statistics));
	cluster.size_distribution.bin = (Statistics *)calloc(psl->nparts, sizeof(Statistics));

	cluster.size_histogram.length=psl->nparts;
	cluster.size_distribution.length=psl->nparts;


	return;
}

void free_all_memory(void){
	// free the memory you allocated in allocate_memory, or other parts in the code

	free(cluster.size_histogram.bin);
	free(cluster.size_distribution.bin);

	free(&copyslice[0]);
	free(&start_slice[0]);
	free(&slice[0]);

	if (sys.sim_type==MC_ALGORITHM){
        free(&cp_slice[0]);
    }

	return;
}
