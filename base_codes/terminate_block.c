#include "path.h"


/*------------------LOCAL FUNCTIONS------------------------------------------*/
void printstatusbmd(Slice *psl);
void count_empty_clusters_monomers(Slice *, int *, int *);
void print_pos_sites(Slice *);
// void printrdf();
// PUT IN VERSION SPECIFIC c file: void print_association_dissociation_times(void);

/*---------------------------------------------------------------------------*/



void terminate_block(Slice *psl) {
    // terminate_block() resides in the outerloop at main.c
    //  a few things are always calculated and printed (to file), irrespective of the specific measurements of the systems you want: 
    /*      1)  the energy (printed to file). In this functions also the information of all bonds is saved
            2)  the current configuration (printed to file)
            3)  the cluster distributions/histogram, see in energy.c how a bond is defined 
            4)  the number of bonds, monomers, and cycli (printed to screen)
        then, depending on the algorithm. Details about it are printed to the screen, e.g. time  in BMD or acceptance ratio's in MC
        Finally, version specific measurements and calculations are performed. Those are listed in version_specific_analysis(psl) in version_specific_functions.c
        
        init: init=1 when you first use terminate_block at the beginning of the program. 
        when init=1, you also save the starting configuration in start_slice (a global variable).
        In some versions of the code, you want to know information about the first snapshot. So this information will always be available in this slice.
    */
        
    static int init=1, teller=0;
    
    // printf("terminate_block\n")
    psl->energy=total_energy(psl);
    printf("         the total energy is %lf\n", psl->energy);
   
    // print to file
    printenergy(psl, ""); 
    conf_output(psl);   

    // printf("identify clusters and count bonds\n"); // identify clusters and count bonds
    int monomer=0,emptycluster=0;
    // always first perform cluster_analysis, then clustersize_identification.
    psl->nclusters = cluster_analysis(psl);
    clustersize_identification(psl);
    count_empty_clusters_monomers(psl, &emptycluster, &monomer);
    int nrings=(psl->nbonds + psl->nclusters)-psl->nparts  ; // quick way to count # cycli . if a ring has formed, the (psl->nbonds + psl->nclusters) > psl->nparts 
    printf("         nbonds= %d       in %d clusters with %d monomers and %d cycli\n", psl->nbonds, psl->nclusters, monomer,nrings);

    switch(sys.sim_type){
        case BMD_ALGORITHM:         
            printstatusbmd(psl); 
            break;
        case MC_ALGORITHM: 
            printstatusmc();
            optimizemc();
            break;
    }
    
    // here the version specific measurements etc
    version_specific_analysis(psl);

    if (init){
        memcpy(start_slice,psl,sizeof(Slice));
        init=0;
    }

    return;
}

void count_empty_clusters_monomers(Slice *psl, int *emptyc, int *nmono){
    // count empty clusters and monoemrs
    int check=0;
    for(int id=0;id<psl->nparts;id++){
        /* add up all clustersizes; it should be equal to npart*/
        check+=cluster.clustersizes[id];
        // 
        if(cluster.clustersizes[id]==0){
            *emptyc +=1;
        }
        if(cluster.clustersizes[id]==1){
            *nmono +=1;
        }
    }
    if (check!=psl->nparts){
        error(" you are counting less particles in clusteranalysis than there are particles ");
    }

    return;

}
    
void print_clusterinformation(Slice *psl ){
    for(int id=0;id<psl->nclusters;id++){
        printf("   ID # %d \\w %d particles:   ", id, cluster.clustersizes[id]);
        for (int n=0;n<cluster.clustersizes[id];n++){
            printf(" %d, ",cluster.pic[id].stack[n]);
        }
        printf(" \n");
    }
} 


void conf_output(Slice *psl) {
    int ipart, isite;
    FILE *fp;


    if((fp=fopen("conf.out","w"))==NULL) {;
        printf("Warning: can not open conf.out\n");
    }
    else {
        fprintf(fp,"%d %.6f  %.6f  %.6f \n", psl->nparts, sys.boxl.x,sys.boxl.y,sys.boxl.z);
        for(ipart=0; ipart<psl->nparts; ipart++) {
            fprintf(fp,"%.6f %.6f %.6f %.6f %.6f %.6f %.6f\n", psl->pts[ipart].r.x, psl->pts[ipart].r.y, psl->pts[ipart].r.z,
                                 psl->pts[ipart].q.q0, psl->pts[ipart].q.q1, psl->pts[ipart].q.q2, psl->pts[ipart].q.q3);
        } }
    fclose(fp);
    return;
}

void finalstat(Slice *psl) {

    int i,istate,jstate,irep,id;
    // printf("starting with finalstat\n");
    printf("\nFinished Simulation\n");


    if(sys.sim_type==0) {
        printstatusbmd(psl); // print the energy

    }
    else if(sys.sim_type==2) {
        printf("\n*** MC stats **\n");
        optimizemc();
        final_printstatusmc_sub(mc_single);

        if(sys.cluster_MC==1){
            final_printstatusmc_sub(mc_cluster);
            final_printstatusmc_sub(mc_mono);
        }
    }
    else{
        error("use only simtype 2 for MC or 0 for bmd");
    }
    
    //  print RDF to file if rdfanalysis is turned on
    // if(sys.rdfanalysis==1){
    //     printrdf();
    // }
    
    return;
}

void print_slice_information(Slice *psl){
    // printing slice information
    int m;
    printf("\n     >>>>> Slice Information <<<<<<\n");
    printf("            energy               %10.10lf\n",psl->energy);
    printf("            nbonds               %d\n",psl->nbonds);
    printf("            nparts               %d\n",psl->nparts);
    if (sys.sim_type==0){
    printf("             c_time               %10.10lf\n",psl->c_time);
    }

    // particle information
    for(int n=0;n<psl->nparts;n++){
        //position 
        printf("\n                 particle      %d     \n",n);
        printf("                 cluster id    %d     \n",psl->pts[n].cluster_id);
        printf("                 ptype         %d     \n",psl->pts[n].ptype);
        printf("                 position    ( %10.5lf  %10.5lf  %10.5lf)     \n",psl->pts[n].r.x,psl->pts[n].r.y,psl->pts[n].r.z);
        printf("                 force       ( %10.5lf  %10.5lf  %10.5lf)     \n",psl->pts[n].f.x,psl->pts[n].f.y,psl->pts[n].f.z);
           
        for(m=0; m<sys.particletype[psl->pts[n].ptype].nsites;m++){
            printf("                 patchvector          ( %10.5lf  %10.5lf  %10.5lf)     \n",psl->pts[n].patchvector[m].x,psl->pts[n].patchvector[m].y,psl->pts[n].patchvector[m].z);
        }
            printf("                 has %d bonds with:           ",psl->pts[n].nbonds);

        for(m=0; m<psl->pts[n].nbonds;m++){
           printf("%d ",psl->pts[n].bonds[m]);
           printf("(EBond=%.15lf),    ",psl->pts[n].bond_energy[m]);
        }
        printf("\n");
    }
    return;
}




void printstatusbmd(Slice *psl) {

    // printf("printstatusbmd\n");
    if(sys.nearest_neighbor) {
        forcecheck_nn(psl);
        energycheck_nn(psl);
    }
    printf("cumulative time  %lf s\n", psl->c_time);

    return;
}



void printenergy(Slice *psl, char ext[]){
    /*Printing the energy to a file. Each printed energy is after ncycle2 steps */
    char value[MAX_VALUE_LENGTH];
    
    switch(sys.sim_type){
        case MC_ALGORITHM:
            snprintf( value, sizeof( value ), "%8.12lf\n", psl->energy );
            break;
        case BMD_ALGORITHM:
            snprintf( value, sizeof( value ), "%8.12lf %8.12lf\n", psl->c_time, psl->energy );
            break;
    }

    write_append_to_file( "energy",  ext, 'a' ,value);

    return;

}


void printing_trajectory(Slice *psl){

    Pts *psi;

    /*saving the coordinates of the X particles in a struct. if bond has broken then print the 1e6 datapoints to a file*/
    
    FILE *coordinatesfile;
    
    double ctime,x,y,z,q0,q1,q2,q3;
    

    if((coordinatesfile=fopen("trajectory.xyz","a"))==NULL) {
            printf("Error with opening \"trajectory.xyz\" file");
            error("stop");
    }
    else{
        ctime=psl->c_time;
        
        for (int ipart=0;ipart<sys.npart;ipart++){
            psi=&psl->pts[ipart];
            x=psi->r.x;
            y=psi->r.y;
            z=psi->r.z;

            q0=psi->q.q0;
            q1=psi->q.q1;
            q2=psi->q.q2;
            q3=psi->q.q3;

            fprintf(coordinatesfile, "%.1f %d %.4lf %.4lf %.3lf %.3lf %.3lf %.3lf %.3lf\n", 
            ctime,ipart,x,y,z,
            q0,q1,q2,q3 );  

        }
        // fprintf(coordinatesfile, "\n");
    }

    fclose(coordinatesfile);

    return;
}


void print_pos_sites(Slice *psl){
    //looop over the particles and print its position and its patches
    // was used for debugging 
    int ipart,p, s;
    Pts *psi;

    for(p=0;p<sys.nparticle_types;p++){
        printf("%d particles have patch: \n",  sys.particletype[p].nparticles);
        for(s=0; s<sys.particletype[p].nsites;s++){
            vprint(sys.particletype[p].site[s]);
        }
    }
    return;
}



void print_StatsLength_to_file(StatsLength *stats_name){
    FILE *file;
    int i, do_print=0;

    // check if there is data in stats_name 
    for(i=0;i<stats_name->length;i++){ 
        if (stats_name->bin[i].n>0){
            do_print=1;
            break;
        }
    }

    if (do_print){
        if ((file = fopen(stats_name->filename,"w"))==NULL){
            printf("%s \n", stats_name->filename);
            error("file can't be opened \n");
        }
        else{
            // printing to file
            for(i=0;i<stats_name->length;i++){ 
                fprintf(file, "%8.12lf %8.12lf %ld\n", stats_name->bin[i].mean, stats_name->bin[i].variance2,stats_name->bin[i].n); 
            }
        }
        fclose(file);
    }

    return;
}



void write_append_to_file(char filename[],  char ext[], char writetype, char value[] ){
    // either (over)writes or appends to file, writetype="w" or "a", resp.
    // extension is maybe extra part of filename you want to have
     FILE *fp;
    char filename_full[205];

    // checks if given strings are not too long.
    if ((strlen(filename) >= MAX_FILENAME_LENGTH) || (strlen(ext) >= MAX_EXT_LENGTH) || (strlen(value) >= MAX_VALUE_LENGTH)){
        error(" In write_append_to_file():  filename, ext, or value exceeds maximum length");
    }

    // the full filename will be composed of 2 parts: filename+ext+".out"
    snprintf(filename_full, sizeof(filename_full), "%s%s.out", filename, ext);

    // printing the value to file, choose writing (w) or appending (a)
    if (writetype != 'w' && writetype != 'a') {
        printf("Error: Invalid write type '%c'. Use 'w' for writing or 'a' for appending.\n", writetype);
        return;
    }

    if ((fp = fopen(filename_full, writetype == 'w' ? "w" : "a")) == NULL) {
        printf("Error: Unable to open file '%s' for %s mode.\n", filename_full, writetype == 'w' ? "writing" : "appending");
        return;
    }

    // print to file
    if (fprintf(fp, "%s\n", value) < 0) {
        printf("Error: Failed to write data to file '%s'.\n", filename_full);
        fclose(fp);
        return;
    }
    // close the file
    fclose(fp); 
    return;
}





