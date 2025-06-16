#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include "path.h"

void printstatusbmd(Slice *psl);
void count_empty_clusters_monomers(Slice *, int *, int *);

// here all the functino for terminate block

void terminate_block(Slice *psl) {

    int count_npart=0, monomer=0,emptycluster=0;
    int n, id,nbonds_old,nbonds_new;
        
    static int init=1;
    // printf("terminate_block\n"); 


    nbonds_old=psl->nbonds;
    psl->energy=total_energy(psl);
    printf("the total energy is %lf\n", psl->energy);
    nbonds_new=psl->nbonds;
    
    printenergy(psl);
    conf_output(psl);

    
    // printf("identify clusters and count bonds\n"); // identify clusters and count bonds
    psl->nclusters = cluster_analysis(psl);
    // printf("identify clustersize_identification\n"); // identify clusters and count bonds
    clustersize_identification(psl);
    // printf("done\n");
    // count_empty_clusters_monomers(psl, &emptycluster, &monomer);

    // int nrings=(psl->nbonds + psl->nclusters)-psl->nparts  ;
    // printf("         nbonds= %d       in %d clusters with %d monomers and %d cycli\n", psl->nbonds, psl->nclusters, monomer,nrings);

    if (analysis.local_density) local_density(psl);
    


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
        // printf(" id %d has %d particles\n",id,psl->cluster.clustersize[id] );
        /* add up all clustersizes; it should be equal to npart*/
        check+=psl->cluster.clustersize[id];
        // 
        if(psl->cluster.clustersize[id]==0){
            *emptyc +=1;
        }
        if(psl->cluster.clustersize[id]==1){
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
        printf("   ID # %d \\w %d particles:   ", id, psl->cluster.clustersize[id]);
        for (int n=0;n<psl->cluster.clustersize[id];n++){
            printf(" %d, ",psl->cluster.pic[id].stack[n]);
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

void print_statistics_file(StatsLength *stats_name, int length){
    FILE *file;
    char *pt,line[NPART], filename[100];
    int i=0, do_print=0;

    for(i=0;i<length;i++){ /*the length histogram starts with length=1; the minimum length of a chain. therefore i+1*/ 
        if (stats_name->bin[i].n>0){
            // printf("at bin %d there are more than zero measurements %ld\n",i,stats_name->bin[i].n);
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
            for(i=0;i<length;i++){ /*the length histogram starts with length=1; the minimum length of a chain. therefore i+1*/ 
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
    // overwrite the file
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

