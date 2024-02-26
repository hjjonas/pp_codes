/*---------------------------------------------------------------------------------------------------------*/
/*      here are analysis function that are NOT specific to a certain simulation or project */
/*---------------------------------------------------------------------------------------------------------*/

#include "path.h"

/*------------------LOCAL FUNCTIONS------------------------------------*/
void DFS(Slice *,int  ,IntArray *, Slice *);
/*---------------------------------------------------------------------------*/



void linking_all_cluster(Slice *psl){
    // printf(" linking_all_cluster \n");

    psl->nclusters = cluster_analysis(psl); 
    clustersize_identification(psl);
    for (int n=0; n<psl->nclusters;n++){
        linked_structure_icluster(psl,n,  cluster.pic[n].stack[0]); 
    }
    return;
    
}

void add_bond_information(Slice *psl,int ipart,int jpart,double Erep,double Ebond,double s){
    // add the information about the bond to ipart, it makes bond with jpart
    // NOTE IT DOES NOT DO THE REVESED , so ipart is not saved to jpart structure
    Pts *psi;
    psi=&psl->pts[ipart];

    psi->bonds[psi->nbonds]=jpart;
    psi->bond_energy[psi->nbonds]= Erep+Ebond;
    psi->bond_distance[psi->nbonds]=s;
    psi->nbonds++; 
    // printf("particle %d (has now %d bonds) makes bond iwth %d. \n",ipart, psi->nbonds,jpart);
    return;
}

void linked_structure_icluster(Slice *psl, int icluster, int ipart_ref){
    // if a structure crossed the pbc, rotating it is not straightforward. 
    // with this fucntion you first check if structure crosses pbc ; are there particles 1 isgma from the border? in the 2D plane
    int i,ipart,jpart;
    int clusteri_size=cluster.clustersizes[icluster];
    vector dr; 
  
    Slice *copyslice = malloc(sizeof(Slice));
    memcpy(copyslice, psl, sizeof(Slice));

    Picl pici=cluster.pic[icluster];
    IntArray particlesleft_list; //

    //initiate particlesleft_list (filled)
    initIntArray(&particlesleft_list,  (size_t) clusteri_size);
    for(i=0; i<clusteri_size; i++){         
        insertIntArray(&particlesleft_list, pici.stack[i]);
    }
    
    // printf("whole array of particles_left\n");
    // you will need to loop over the bonds, because the boundaries are crossed by the structure
    // start from ipart_ref; 
    ipart=ipart_ref; 
    removeElementXIntArray( &particlesleft_list ,  ipart); //remove current(=ipart) from particlesleft list

    // walk untill you reach particle with 1 bond or run into a particles thats not in the list. 
    for ( i=0; i<psl->pts[ipart].nbonds ; i++ ){
        jpart=psl->pts[ipart].bonds[i];

        // printf(" perform DFS with   start from particle %d \n", start);
        if (checkElementXIntArray(&particlesleft_list, jpart )==1){
            // new position next
            dr=particles_vector(psl,  ipart, jpart); //return vector from ipart to jpart
            vector_add(copyslice->pts[ipart].r,dr,copyslice->pts[jpart].r); // of course dont do pbc after this step! 

            DFS( psl,jpart,  &particlesleft_list ,copyslice );
        }
    }
    
    // when done, free the memory!
    freeIntArray(&particlesleft_list);

    //copy new structure to psl;
    memcpy(psl,copyslice,  sizeof(Slice));
    free(copyslice);
    return;
}

void DFS(Slice *psl,int ipart, IntArray *particlesleft_list, Slice *copyslice ){
    //  Depth First Search  adapted from https://www.codewithharry.com/videos/data-structures-and-algorithms-in-hindi-89/
    //  By HJ Jonas oktober 2022
    //  it walks over the particles and performs DPF
    //  it makes smart use of a nested loop, by recalling DFS inside DFS, to walk over the structure 
    //  see section 2.5.3 DepthFirstSearch of thesis HJJ
    vector dr;
    int jpart,j_inlist;

    // remove ipart from list, else you will find the bond back to ipart. Now you will walk forward to next bond
    removeElementXIntArray(particlesleft_list,  ipart );

    // loop over the bound particles of ipart, to see if you have visited them
    for (int n = 0; n < psl->pts[ipart].nbonds; n++){   
        jpart=psl->pts[ipart].bonds[n]; // the bound particles to ipart;

        //  did you already visit jpart?
        j_inlist=checkElementXIntArray(particlesleft_list, jpart );
        if(j_inlist==1){
            // new position jpart, and perform DFS on jpart
            dr=particles_vector(psl,  ipart, jpart ); // points from i to j
            vector_add(copyslice->pts[ipart].r,dr,copyslice->pts[jpart].r);     
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
    dist=dr - (sys.particletype[psl->pts[i].ptype].radius+sys.particletype[psl->pts[j].ptype].radius);
    
    //return
    *distance=dist;
    *length=dr;
    *vec=bond_vec;
       
    
    return;
}

double particles_distance(Slice *psl, int i, int j){
    /* calculates and returns the surface-surface distance between two particles. specify with pbc if pbc should be on or off*/

    double dr2, r;
    vector dr;

    dr = particles_vector(psl, i, j);
    dr2 =vector_inp(dr,dr);
    r=sqrt(dr2) - (sys.particletype[psl->pts[i].ptype].radius+sys.particletype[psl->pts[j].ptype].radius);
    
    return r;
}

vector particles_vector(Slice *psl, int i, int j){
    /* calculates and returns the vector rij(=ri-rj) between two particles. 
        return vector from i pointing to j*/
    vector dr;

    vector_minus(psl->pts[j].r,psl->pts[i].r,dr); 
    pbc(dr, sys.boxl); 
    
    return dr;  
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
    return ;              
}

int cluster_analysis(Slice *psl){
    /* color all particles white, loop over all other particles (j) color them grey; find bonding particles (with dr and SiSj)*/
    /*particles are in a cluster if they have a bond energy <sys.bond_cutoffE*/
    // code by PGB
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

    nwhite = psl->nparts;
    for(i=0;i<psl->nparts;i++) {
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
                        /*            printf("deeltje %d wordt zwart,nwhite= %d, ngrey= %d nblack= %d, ntot=%d\n",j,nwhite,ngrey,nblack,nwhite+nblack+ngrey);*/
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

    return ncluster;
}

void clustersize_identification(Slice *psl){
    /* puts all particles in the correct stack if its clusterid = cluster.pic[id].stack
    and counts how many clusters of a certain size there are in cluster.clustersize.*/
    int id,ipart, nclusters=psl->nclusters, n;
    Picl pici;
    int  clustersizes[NPART]={0};

    // printf("loop over all particles en deel ze in in een picl\n" );
    //loop over all particles en deel ze in in een picl 
    for(ipart=0;ipart<psl->nparts;ipart++){
        id = psl->pts[ipart].cluster_id;
        cluster.pic[id].stack[clustersizes[id]]=ipart;
        clustersizes[id]++; // this is a list which tells you how many particles the cluster with clusterid "id" is 
    }
    memcpy(cluster.clustersizes, clustersizes, sizeof(clustersizes));

    return;
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
