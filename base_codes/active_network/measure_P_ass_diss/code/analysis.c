#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "path.h"

void DFS(Slice *,int  ,IntArray *, Slice *);

void add_bond_information(Slice *psl,int ipart, int isite, int jpart, double Erep,double Ebond,double s){
    // add the information about the bond to ipart, it makes bond with jpart
    // NOTE IT DOES NOT DO THE REVESED , so ipart is not saved to jpart structure
    Pts *psi;
    psi=&psl->pts[ipart];

    psi->bonds[psi->nbonds]=jpart;
    psi->bound_site[psi->nbonds]=isite;
    psi->bond_energy[psi->nbonds]= Erep+Ebond;
    psi->bond_distance[psi->nbonds]=s;
    psi->nbonds++; 
    // printf("particle %d (has now %d bonds) makes bond iwth %d. \n",ipart, psi->nbonds,jpart);
    return;

}

void linked_structure_icluster(Slice *psl, int icluster, int ipart_ref){
    // if a structure crossed tha pbc, rotating it is straightforward. 
    // with this fucntion you first check if structure crosses pbc ; are there particles 1 isgma from the border? in the 2D plane
    int i,b,current,ipart,boundary_crossing=0,nbonds;
    int particlesleft=psl->cluster.clustersize[icluster];
    vector loc_prev;
    int clusteri_size=psl->cluster.clustersize[icluster];
    vector dr;
    Pts newpos[NPART],p_ref;
    

    // printf(" copy psl to newpos\n");
    memcpy(newpos, psl->pts, sizeof(psl->pts));
    p_ref = psl->pts[ipart_ref];

    // vprint(p_ref.r);
    //move the structure such that p_ref is in the center, while checking if the new structure is near the boundary after pbc
    int edge=0;
    for(i=0; i<clusteri_size; i++){
        ipart= psl->cluster.pic[icluster].stack[i];
        
        pbc(newpos[ipart].r,sys.boxl);
        
        if ((0.5*sys.boxl.x-fabs(newpos[ipart].r.x)<1.5 )| (0.5*sys.boxl.y-fabs(newpos[ipart].r.y)<1.5)| (0.5*sys.boxl.z-fabs(newpos[ipart].r.z)<1.5) ){
            edge=1;
            // vprint(newpos[ipart].r);
            // vector_add(newpos[ipart].r,p_ref.r,newpos[ipart].r); // so newpos[ipart].r is the chain with p_ref in (0,0,0)
            // printf(" part of structure (particle %d ) near the boundary\n",ipart);
            break; // the structure is near the boundary and may cross it
        }
    }

    if (edge==0){
        // printf(" no particles at edge, newpos is correct as it is\n");
        //no particles at edge, newpos is correct as it is
        memcpy( psl->pts, newpos, sizeof(psl->pts));
        return;
    }

    Slice *copyslice = malloc(sizeof(Slice));
    memcpy(copyslice, psl, sizeof(Slice));

    if(edge==1){
        // printf(" edge =1 \n");
        Picl pici=psl->cluster.pic[icluster];
        IntArray particlesleft_list; //,branchlist; 
        int current;
        vector dr;

        //initiate particlesleft_list (filled)
        initIntArray(&particlesleft_list,  (size_t) clusteri_size);
        for(i=0; i<clusteri_size; i++){         
            insertIntArray(&particlesleft_list, pici.stack[i]);
        }
        
        // printf("whole array of particles_left\n");
        // you will need to loop over the bonds, because the boundaries are crossed by the structure
        //start from ipart_ref; 
        current=ipart_ref; 
        // newpos[current].r=psl->pts[current].r; // the starting point/particle, check if it is a " branchpoint", you need to go backwards too
        removeElementXIntArray( &particlesleft_list ,  current); //remove current from particlesleft

        // walk untill you reach particle with 1 bond or run into a particles thats not in the list. 
        for ( i=0; i<psl->pts[current].nbonds ; i++ ){
            int next=psl->pts[current].bonds[i];

            // printf(" perform DFS with   start from particle %d \n", start);
            if (checkElementXIntArray(&particlesleft_list, next )==1){
                // new position next
                dr=particles_vector(psl, current, next); //check if order of( current, next) is correct (direction of vector)
                vector_minus(copyslice->pts[current].r,dr,copyslice->pts[next].r); // of course dont do pbc after this step! 

                DFS( psl,next,  &particlesleft_list ,copyslice );
            }
        }
        
        // when done, free the memory!!
        freeIntArray(&particlesleft_list);
    }

    //copy new structure to psl;
    memcpy(psl,copyslice,  sizeof(Slice));
    free(copyslice);
    return;
}

void DFS(Slice *psl,int ipart, IntArray *particlesleft_list, Slice *copyslice ){
    //  Depth First Search  adapted from https://www.codewithharry.com/videos/data-structures-and-algorithms-in-hindi-89/
    //  By HJ Jonas oktober 2022
    //  it walks over the particles and performs DPF
    //  it makes smart use of a nested loop, by recalling DFS inside DFS

    vector dr;
    // remove ipart from list, else you will find the bond back to ipart. Now you will walk forward to next bond
    removeElementXIntArray(particlesleft_list,  ipart );

    // printf(" particle %d has %d bonds:  \n", ipart, psl->pts[ipart].nbonds);
    // loop over the bound particles of ipart, to see if you have visited them
    for (int n = 0; n < psl->pts[ipart].nbonds; n++){   
        int jpart=psl->pts[ipart].bonds[n]; // the bound particles to ipart;

        // do ipart and jpart make a bond? and is jpart still in the particlesleft_list list
        int make_bond=bond_check(psl,ipart,jpart);
        int j_inlist=checkElementXIntArray(particlesleft_list, jpart );

        // printf(" particle %d and %d make bond (make_bond=%d), && (j_inlist=%d ==1?).  \n", ipart,jpart,make_bond,j_inlist);
        if((make_bond==1) && (j_inlist==1)){
            // printf("    particle %d and %d make bond, and jpart=%d was not visited yet.  \n", ipart,jpart,jpart);

            // new position jpart
            dr=particles_vector(psl, ipart, jpart); //check if order of( current, next) is correct (direction of vector)
            vector_minus(copyslice->pts[ipart].r,dr,copyslice->pts[jpart].r);
            
            DFS( psl,jpart,  particlesleft_list,  copyslice );
        }
    }

    return;
}



int bond_check(Slice *psl, int i, int j){
    /* check whether there exists a bond between particle i and j. Based on the treshold value sys.bond_cutoffE
    returns 1 for yes, returns 0 for no*/
    Pts *psi,*psj; 
    int n;  

    psi=&psl->pts[i];
    psj=&psl->pts[j];    

    for (n=0;n<psi->nbonds;n++){
        if (psi->bonds[n]==j){
            return 1;
        }
    }
    for (n=0;n<psj->nbonds;n++){
        if (psj->bonds[n]==i){
            return 1;
        }
    }

    return 0;
}

void particles_distance_length_vector(Slice *psl, int i, int j, double *distance, double *length,vector *vec){
    // returns :
    //   the surface, surface distance beweteen particles (dist),  
    //   center-to-center vector (vec), and 
    //   length of  center-to-center vector (length)

    double dr2, r;
    double dr, dist;
    vector bond_vec;


    bond_vec = particles_vector(psl, i, j);
    dr2 =vector_inp(bond_vec,bond_vec);
    dr=sqrt(dr2);
    dist=dr - sys.particletype[psl->pts[i].ptype].radius-sys.particletype[psl->pts[j].ptype].radius;
    
    //return
    *distance=dist;
    *length=dr;
    *vec=bond_vec;
       
    
    return;
}

double particles_distance(Slice *psl, int i, int j){
    /* calculates and returns the distance between two particles. specify with pbc if pbc should be on or off*/

    double dr2, r;
    vector dr;

    dr = particles_vector(psl, i, j);
    dr2 =vector_inp(dr,dr);
    r=sqrt(dr2) - sys.particletype[psl->pts[i].ptype].radius-sys.particletype[psl->pts[j].ptype].radius;
    
    return r;
}



vector particles_vector(Slice *psl, int i, int j){
    /* calculates and returns the vector rij(=ri-rj) between two particles. specify with pbc if pbc should be on(pbc=1) or off(pbc=0)*/
    vector dr;

    vector_minus(psl->pts[i].r,psl->pts[j].r,dr); /* dr of particle j and k */
    pbc(dr, sys.boxl); 
    
    return dr;  
}

void clustersize_freq_update(Slice *psl){
    /*this function saves the frequency of occurence of the clusterlengths  */

    int i,l;
    int   freq_clusterlength[NPART]={0};
    
    for(i=0;i<psl->nclusters;i++){
        freq_clusterlength[psl->cluster.clustersize[i]-1]++; // use -1 here because size one is the smallest size
    } 

    for(l=0;l<psl->nparts;l++){
        /*loops over the lengths which has maximumvalue of psl->nparts*/
        // update_average(chain.length_histogram[i], freq_clusterlength[i]);
        
        running_statistics(&psl->cluster.size_histogram.bin[l], freq_clusterlength[l]);
        running_statistics(&psl->cluster.size_distribution.bin[l], (double)freq_clusterlength[l]/psl->nclusters);
    } 

    return; 
}



void check_maxbonds(Slice *psl){
    /*checks if any particle exeeds the maximum number of bonds (more than one per site)*/
    int ipart;
    Pts *psi;

    for(ipart=0;ipart<psl->nparts;ipart++){
        psi=&psl->pts[ipart];
        if (sys.particletype[psi->ptype].nsites==NSITES){
            //if sys.particletype[psi->ptype].nsites==NSITES; the particle is isotropic and may have more than NSITES bonds
            continue;
        }
        else if((psi->nbonds>sys.particletype[psi->ptype].nsites)){
            dprint(ipart);
            dprint(psi->nbonds);
            dprint(sys.particletype[psi->ptype].nsites);
            error("too many bonds for the particle");
        }
    }
    return;              
}

int cluster_analysis(Slice *psl){
    /* color all particles white, loop over all other particles (j) color them grey; find bonding particles (with dr and SiSj)*/
    /*particles are in a cluster if they have a bond energy <sys.bond_cutoffE*/
    Pts *psj;
    vector dr,norm;
    int label[psl->nparts],nwhite,notallgreysgone,ngrey,ncluster,nblack,i,j,k, WHITE=1, BLACK=2, GREY=3, id;
    double r,dr2,s,r_radii,  potential_energy_P_d=0;
   

    //  count=0;
    notallgreysgone=0;
    ngrey=nblack=ncluster = 0;
    for(i=0;i<psl->nparts;i++) {
        label[i] = WHITE;
    }

    // printf(" >> inside cluster_analysis\n");

    nwhite = psl->nparts;

    for(i=0;i<psl->nparts;i++) {
        // dprint(i);
        if (label[i] == WHITE) {
            nwhite--;
            ngrey =1;
            label[i]= GREY;
            do {
                notallgreysgone = 0;
                for (j=0;j<psl->nparts;j++) {
                    if (label[j]==GREY) {
                        notallgreysgone = 1;
                        ngrey--;
                        nblack++;
                        label[j]=BLACK;
                        psj = &psl->pts[j];
                        psj->cluster_id=ncluster; /* give particle j the cluster_id of particle i*/
                                    // printf("deeltje %d wordt zwart,nwhite= %d, ngrey= %d nblack= %d, ntot=%d\n",j,nwhite,ngrey,nblack,nwhite+nblack+ngrey);
                        for (k=i;k<psl->nparts;k++){
                            if ((k!=j) && (label[k] == WHITE)) {
                                // for bond_check to work, you should run total_energy before the cluster analyssis! 
                                if(bond_check(psl, j,k)){
                                    label[k]=GREY;
                                    nwhite--;
                                    ngrey++;
                                }
                            }
                        }
                    }
                }
                if (ngrey==0)   notallgreysgone=0;
            } while (notallgreysgone);
            ncluster++;
            if (nwhite==0) i= psl->nparts;
        }
    }

    // printf("the end. there are %d clusters \n",ncluster );

    return ncluster;
}

void clustersize_identification(Slice *psl){
    /* puts all particles in the correct stack if its clusterid = psl->cluster.pic[id].stack
    and counts how many clusters of a certain size there are in psl->cluster.clustersize.*/
    int id,ipart, nclusters=psl->nclusters, n;
    Picl pici;
    int  clustersize[NPART]={0};

    // printf("loop over all particles en deel ze in in een picl\n" );
    //loop over all particles en deel ze in in een picl 
    for(ipart=0;ipart<psl->nparts;ipart++){
        id = psl->pts[ipart].cluster_id;
        psl->cluster.pic[id].stack[clustersize[id]]=ipart;
        clustersize[id]++;
    }
    memcpy(psl->cluster.clustersize, clustersize, sizeof(clustersize));

    return;
}


double Rg2_radius_of_gyration(Slice *psl, int icluster){
    /*https://en.wikipedia.org/wiki/Radius_of_gyration
    : Mathematically the radius of gyration is the root mean square distance of the object's parts from either 
    its center of mass or a given axis, depending on the relevant application.

    Rg^2 = 1/N \sum_k=1^N( \vec{r}_k - \vec{r}_{mean})^2
    where N is number of particles (in cluster)
    vec r_k is the position of particle k
    vec r_mean is the average postiion of the particles in the cluster
    */
    int ipart,ipart_ref, clustersize,n;
    vector ipart_r,rdif;
    double Rg2=0,r2;

    ipart_ref= psl->cluster.pic[icluster].stack[0];
    clustersize=psl->cluster.clustersize[icluster];

    linked_structure_icluster(psl,icluster,  ipart_ref);
    ipart_r=psl->pts[ipart_ref].r;

    for(n=0;n<clustersize;n++){
        ipart=psl->cluster.pic[icluster].stack[n];
        vector_minus(psl->pts[ipart].r,ipart_r,rdif);
        r2=vector_inp(rdif,rdif);
        Rg2+=r2;
    }

    return (Rg2/clustersize);

}


int particle_in_wall(Slice *psl, int ipart){
    /* returns 1 if particle is in the wall, if not, it returns 0. The wall is defined at 1*/
    if(psl->pts[ipart].r.z<=sys.particletype[psl->pts[ipart].ptype].radius){
        return 1;
    }
    else{
        return 0;
    }
}

void reset_running_statistics(Statistics *stats){
    /* reset the statistics; set values to zero  */

    stats->mean=0.;
    stats->variance2=0.;
    stats->n=0;
    return;
}

double running_variance2(double x_new, double s_old2, double u_old, double u_new, long n_old){
    /*  x_new is the new measurement,
        s_old2 is the old variance^2 and stored in chain.e2e2_sigma2
        u_old is the old mean
        u_new is the new mean
        s_old2 is stored in chain.e2e2_sigma2[chainlength];
        (N-1)S_new^2 =  (N-2)S_old^2 +(x_new-u_new)(x_new-u_old)
        
        return  S_{new}^2 */

    double s_new2, s_new2_N_1;
    long n_new=n_old+1;

    s_new2_N_1 = s_old2*n_old + (x_new-u_new)*(x_new-u_old);
    s_new2 = s_new2_N_1 /(n_new);
    return s_new2;
}

double running_mean(double u_old, double x_new, long n_old){
    /*  x_new is the new measurement,
        u_old is the old mean
        n_old is the old number of measurements;
        n_new = n_old+1
        u_new = u_old + (x_new-u_old)/n_new
        
        return  n_new */
    long n_new = n_old+1;
    double u_new;

    u_new = u_old + (x_new-u_old)/n_new;
    return u_new;
}

void running_statistics(Statistics *oldstats, double x_new){
    /* updates the statistics
    it needs the old statistics (oldstats) and the value of the new measurement
    the new mean value and the new variance is calculated and the statistics are updated.  */

    double u_old =oldstats->mean, s_old2 =oldstats->variance2;
    long n_old = oldstats->n;
    double u_new, s_new2;

    u_new = running_mean(u_old,  x_new,  n_old);
    s_new2 = running_variance2( x_new,  s_old2,  u_old,  u_new,  n_old);

    oldstats->mean=u_new;
    oldstats->variance2=s_new2;
    oldstats->n++;
}
