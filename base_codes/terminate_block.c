#include "path.h"


/*------------------LOCAL FUNCTIONS------------------------------------------*/
void printstatusbmd(Slice *psl);
void count_empty_clusters_monomers(Slice *, int *, int *);
void print_pos_sites(Slice *);
void print_energy_s_theta1_theta2(Slice *);
void print_statistics_N1(Statistics , char [100]);
void print_statistics_file(StatsLength *,Slice *);
// void printrdf();
// PUT IN VERSION SPECIFIC c file: void print_association_dissociation_times(void);

/*---------------------------------------------------------------------------*/



void terminate_block(Slice *psl) {

    int monomer=0,emptycluster=0;
    int nbonds_old,nbonds_new;
        
    static int init=1, teller=0;
    // printf("terminate_block\n"); 


    nbonds_old=psl->nbonds;
    psl->energy=total_energy(psl);
    printf("the total energy is %lf\n", psl->energy);
    nbonds_new=psl->nbonds;
    
    printenergy(psl);
    conf_output(psl);   
    
    // If printing trajectory is specified in parameters
    if ((init==0) && (analysis.print_trajectory))  printing_trajectory(psl); 


    // printf("identify clusters and count bonds\n"); // identify clusters and count bonds
    psl->nclusters = cluster_analysis(psl);
    clustersize_identification(psl);
    count_empty_clusters_monomers(psl, &emptycluster, &monomer);

    int nrings=(psl->nbonds + psl->nclusters)-psl->nparts  ;
    printf("         nbonds= %d       in %d clusters with %d monomers and %d cycli\n", psl->nbonds, psl->nclusters, monomer,nrings);

    switch(sys.sim_type){
        case BMD_ALGORITHM:         
            printstatusbmd(psl); // print the energy
            break;
        case MC_ALGORITHM: 
            printstatusmc();
            optimizemc();
            break;
    }



    if((cluster.analysis==1)  && (init==0)){
        clustersize_freq_update(psl);
        // print_statistics_file(&cluster.size_histogram,psl);
        // print_statistics_file(&cluster.size_distribution,psl);
        // local_density(psl);
        // N_TPP_bonds(psl);
        // bond_probability(psl);
    }

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
        // printf(" id %d has %d particles\n",id,cluster.clustersize[id] );
        /* add up all clustersizes; it should be equal to npart*/
        check+=cluster.clustersize[id];
        // 
        if(cluster.clustersize[id]==0){
            *emptyc +=1;
        }
        if(cluster.clustersize[id]==1){
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
        printf("   ID # %d \\w %d particles:   ", id, cluster.clustersize[id]);
        for (int n=0;n<cluster.clustersize[id];n++){
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

void print_adjacency_matrix_coordinates(Slice *psl){
    /* the adjacency matrix is a NxN matrix with zero's and one's
    0 if no bond between particle i and j, and a 1 if there is a bond.
    If the total energy is calculated, the bonds are tracked.  see potential_energy(Slice *psl)
    We will first try to print the whole complete adjacency matrix.
    And evaluate later if the files are not getting too big. THere is a lot of redundant data
    We might want to only print which particles have a bond, and later construct the adjacency matrix in python
    */
    int ipart, jpart,n;
    FILE *fp;
    Pts *psi;

    if((fp=fopen("adjacency_matrix_coordinates.out","a"))==NULL) {;
        printf("Warning: can not open adjacency_matrix.out\n");
    }
    else {
        /*we loop over the particles. 
        in psi->bonds[psi->nbonds]=jpart the particle numbers are printed which have a bond with psi
        start with n=0, there are no bond found with ipart yet*/
        for(ipart=0; ipart<psl->nparts; ipart++) {
            psi = &psl->pts[ipart];
            for(n=0; n<psi->nbonds; n++){
                jpart=psi->bonds[n];
                if(jpart>ipart){
                    fprintf(fp,"%d,%d\n", ipart,jpart); 
                }
            }
        }
    }
    fprintf(fp,"\n");    
    fclose(fp);
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


void printenergy_warmup(Slice *psl){
    /*Printing the energy to a file. Each printed energy is after ncycle2 steps */
    FILE *energyfile;
    char energyfilename[100];
    char* a = "energy_warmup";
    char* extension = ".out";


    
    snprintf( energyfilename, sizeof( energyfilename ), "%s%s", a,  extension );
    energyfile = fopen(energyfilename,"a");


    if (energyfile == NULL){
            printf("Error with opening \"energy_warmup.out\" file");
    }
    else{    
        if(sys.sim_type==0){
            fprintf(energyfile, "%8.12lf %8.12lf\n", psl->c_time, slice[0].energy);  }
        else{
                fprintf(energyfile, "%8.12lf\n", slice[0].energy);   }  
    }
    fclose(energyfile);

    return;

}

void printenergy(Slice *psl){
    /*Printing the energy to a file. Each printed energy is after ncycle2 steps */
    FILE *energyfile;
    char energyfilename[100];
    char* a = "energy";
    char* extension = ".out";

//    printf(" >>>>>>>>>>>>>>>>>>>>PRINGITN ENERGY.OUT\n");
    snprintf( energyfilename, sizeof( energyfilename ), "%s%s", a,  extension );
    energyfile = fopen(energyfilename,"a");


    if (energyfile == NULL){
            printf("Error with opening %s file \n",energyfilename);
            error("Error");
    }
    else{    
        if(sys.sim_type==0){
            fprintf(energyfile, "%8.12lf %8.12lf\n", psl->c_time, slice[0].energy);  }
        else{
                fprintf(energyfile, "%8.12lf\n", slice[0].energy);   }  
    }
    fclose(energyfile);

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

            fprintf(coordinatesfile, "%.1f %d %.4lf %.4lf %.3lf %.3lf %.3lf %.3lf %.3lf\n", ctime,ipart,x,y,z,q0,q1,q2,q3 );  

        }
        // fprintf(coordinatesfile, "\n");
    }

    fclose(coordinatesfile);

    return;
}




void print_energy_s_theta1_theta2(Slice *psl){
    /*Printing the energy to a file. Each printed energy is after ncycle2 steps */
    FILE *energyfile;
    double s, r, cositheta,cosjtheta,cosijtheta, s_treshold=0.5,anglei, anglej;
    int i,j;
    vector rij,u1;
    double Erep, Ec, S;
    


    if((energyfile=fopen("ctime_s_theta1_theta2_energy.out","a"))==NULL) {;
            printf("Error with opening \"print_energy_s_theta1_theta2.out\" file");
    }
    else{  
        for (i=0;i<psl->nparts; i++){
            for (j=i+1;j<psl->nparts; j++){
                particles_distance_length_vector(psl,i,j,&s,&r,&rij); 

                if(s<s_treshold && s>0.0) { 
                    scalar_divide(rij,r,u1);
                    orientation_parameters( psl,  i,  j, u1, &cositheta, &cosjtheta, &cosijtheta);
                    Pts *psi=&psl->pts[i],*psj=&psl->pts[j];
                    
                    anglei=cosangle_to_angle(cositheta);
                    anglej=cosangle_to_angle(cosjtheta);

                    // the energy components of the bond
                    Erep=potential_repulsive_energy_sdist(s);
                    Ec = potential_attractive_energy_sdist(s);
                    S=S_value(cositheta,cosjtheta,0,psi->ptype, psj->ptype);

                    fprintf(energyfile, "%8.6lf %8.12lf %8.12lf %8.12lf %8.12lf %8.12lf %8.12lf %8.12lf\n", psl->c_time,  s, anglei, anglej, Erep+Ec*S, S, Erep, Ec );     
                     
                }
            }
        }
    }
    fclose(energyfile);

    return;

}


void print_pos_sites(Slice *psl){
    //looop over the particles and print its position and its patches
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


// analysis
// void printrdf() {
//     /* this function prints the RDF, which is stored in rdf_averag, to the file "rdf.out". The rdf per slice is calculated in calc_rdf_localCN() */
//     int i;
//     FILE *rdffile;
//     double maxl;
//     MIN(sys.boxl.x,sys.boxl.y,maxl);

//     rdffile = fopen("rdf.out","w");
//     printf("\nPrinting the RDF to rdf.out\n");

//     if (rdffile == NULL){
//             printf("Error with opening \"rdf.out\" file");
//     }
//     else{
//         for(i=0; i<RDFBINS/2; i++){
//             fprintf(rdffile, "%2.4lf  %.5lf\n", i*maxl/RDFBINS, rdf_average[i].sum/rdf_average[i].n);
//         }
//     }
//     fclose(rdffile);
//     return;
// }

void print_statistics_N1(Statistics stats_name, char filename1[100]){
    FILE *file;
    char *pt,line[NPART], filename[100];
    int i=0, dummy;
    char* extension = ".out";
  
    snprintf(filename, sizeof(filename), "%s%s", filename1,extension);
    
    
    if ((file = fopen(filename,"w"))==NULL){
        printf("%s \n", filename);
        error("input: can't be opened \n");
    }
    else{
        // printf("printing statistics to %s\n", filename );
        fprintf(file, "%8.12lf %8.12lf %ld\n",  stats_name.mean, stats_name.variance2,stats_name.n); 
            
    }
    fclose(file);

    return;
}

void print_statistics_file(StatsLength *stats_name, Slice *psl){
    FILE *file;
    char *pt,line[NPART], filename[100];
    int i, do_print=0;

    for(i=0;i<psl->nparts;i++){ /*the length histogram starts with length=1; the minimum length of a chain. therefore i+1*/ 
        if (stats_name->bin[i].n>0){
            // printf("at bin %d there are more than zero measurements %s\n",i,stats_name->bin[i].n)
            do_print=1;
            break;
        }
    }

    if (do_print){
        memcpy(filename,stats_name->filename,sizeof(filename));  
        
        if ((file = fopen(filename,"w"))==NULL){
            printf("%s \n", filename);
            error("input: can't be opened \n");
        }
        else{
            // printf("printing statistics to %s\n", filename );
            for(i=0;i<psl->nparts;i++){ /*the length histogram starts with length=1; the minimum length of a chain. therefore i+1*/ 
                fprintf(file, "%d %8.12lf %8.12lf %ld\n", i+1, stats_name->bin[i].mean, stats_name->bin[i].variance2,stats_name->bin[i].n); 
                // printf("%d %8.12lf %8.12lf %d\n", i, stats_name->length[i].mean, stats_name->length[i].variance2,stats_name->length[i].n);
            }
        }
        fclose(file);
    }

    return;
}

void append_to_file(char filename[100], int ipart, double value){
    FILE *fp;
    char filename_full[150];
    char* extension = ".out";
    
    // the full filname will be composed of 2 parts: filename+"_bondX.out"
    snprintf(   filename_full, sizeof( filename_full ), "%s_ptype%d%s", filename,ipart, extension );

    if((fp=fopen(filename_full,"a"))==NULL) {; // 
        printf("Warning: not able to append to %s\n", filename_full);
        error("stop");
    }
    else { fprintf(fp,"%.10lf\n",value); }
    fclose(fp);

    return;
}

void print_to_file(char filename[100], int ipart, double value){
    FILE *fp;
    char filename_full[150];
    char* extension = ".out";
    
    // the full filname will be composed of 2 parts: filename+"_bondX.out"
    snprintf(   filename_full, sizeof( filename_full ), "%s_ptype%d%s", filename,ipart, extension );


    // printing the S values to file
    if((fp=fopen(filename_full,"w"))==NULL) {; // 
        printf("Warning: not able to write to %s\n", filename_full);
        error("stop");
    }
    else { fprintf(fp,"%.10lf\n",value); }
    fclose(fp);

    return;
}

void append_to_file2(char filename[100], char value[1000], char extension[100]){
    // append the file, no overwriting
    // extension is maybe extra part of filename you want to have
    FILE *fp;

    char filename_full[250];
    
    
    // the full filname will be composed of 2 parts: filename+extension+".out"
    snprintf(   filename_full, sizeof( filename_full ), "%s%s.out", filename, extension );


    if((fp=fopen(filename_full,"a"))==NULL) {; // 
        printf("Warning: not able to append to %s\n", filename_full);
        error("stop");
    }
    else { fprintf(fp,"%s\n",value); }
    fclose(fp);

    return;
}

void print_to_file2(char filename[100], char value[1000], char extension[100]){
    // overwrite the file
    // extension is maybe extra part of filename you want to have
    FILE *fp;

    char filename_full[250];
    
    
    // the full filname will be composed of 2 parts: filename+extension+".out"
    snprintf(   filename_full, sizeof( filename_full ), "%s%s.out", filename, extension );


    // printing the S values to file
    if((fp=fopen(filename_full,"w"))==NULL) {; // 
        printf("Warning: not able to write to %s\n", filename_full);
        error("stop");
    }
    else { fprintf(fp,"%s\n",value); }
    fclose(fp);

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



