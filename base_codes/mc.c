#include "path.h"
/* ___________________ THE MONTE CARLO ALGORITHM __________________________*/

/*------------------LOCAL FUNCTIONS------------------------------------------*/
void copy_clusterparticles( Slice * ,Slice *,int );
int check_internal_energy(Slice *,Slice *, int , char [100]);
int rotate_monocluster(Slice *, int );
int energy_divergence_check(Slice *, char [50]);
int energy_divergence_check_ipart(Slice *, Pts, int, char [50]);
void single_particle_move(Slice *, int );
int rotatepart_cluster_Ekparts(Slice *);
int translatepart_cluster_Ekparts(Slice *);
void propagate_mc( Slice * );
void cluster_propagate_mc(Slice * );
void printstatusmc_sub(MC *);
void tailflipping(Slice *);
void DFS_tailflip(Slice *,int , IntArray *, Slice * , tensor , quaternion );
/*---------------------------------------------------------------------------*/


void mccycle(Slice *psl) {
    /*sys.sim_type==MC_ALGORITHM
    This function, mccycle, is a Monte Carlo cycle which performs a series of 
    single particle moves or cluster moves depending on a random number generator. 
    The function takes a pointer to a Slice struct as its only argument.*/

    int i,ipart,istep, nclusters,nbonds_before ;
    double which, dr2,r; 
    int nnlupdate=0,j,k,Edivcheck;
    vector dr;
 
    
    /* first identify the clusters; the number of clusters will not change in a clustermove, only the positions of them*/
    /* we first make the clusters with cluster_analysis(psl), followed by identifying the lengths (with clustersize_identification), and linked-list 
    such that the clusters are available via head[chainid] and they are sequenced */

    // printf("nbonds_after %d \n", nbonds_after);
    // printf("start MC\n");
    if(sys.cluster_MC==1){
    // printf("start cluster_MC\n");

        
        for(istep=0;istep<10;istep++){

            /* choose single or cluster movement*/
            which = RandomNumber();

            /* cluster moves, where nbonds and nclusters is not allowed not change*/
            /* the number of cluster should stay constant during clustermoves due to detailed balance*/
            if(which<.05){
                total_energy(psl); // it's necessary to do total_energy before knowing nbonds 
                // printf(" clustermoves, starting iwth %d bonds \n ",psl->nbonds);
                nbonds_before= psl->nbonds;
                // copy cluster if things went wrong to go back in time.
                memcpy(copyslice,psl,sizeof(Slice));

                /* expensive function, do this function only when the number of bonds has changed. */
                if (cluster.update==1){
                    // printf("****performing cluster analysis and linking*****\n");
                    linking_all_cluster(psl);
                    cluster.update=0;

                    // printf("  there are %d clusters \n",psl->nclusters );
                }
                
                /* *** the cluster move*** */
                cluster_propagate_mc(psl);

                // printf("a cluster move is performed, the energy is         %lf\n", psl->energy);
                Edivcheck=energy_divergence_check(psl,"cluster move");
                if (Edivcheck>-1){
                    error("something went wrong in cluster moves, to back to snapshot before cluster moves\n");

                    memcpy(psl,copyslice,sizeof(Slice));
                    psl->energy=total_energy(psl);
                    psl->nclusters = cluster_analysis(psl); 
                    clustersize_identification(psl);
                    nbonds_before= psl->nbonds;
                    linking_all_cluster(psl);
                } //else if (Edivcheck==-1){ everything is ok }

                /*check nclusters is constant, it is constant if nbonds is constant*/
                total_energy(psl);

                // printf(" clustermoves, ending iwth %d bonds \n ",psl->nbonds);
                if(nbonds_before != psl->nbonds){
                    int totbonds=0;
                    for(ipart=0;ipart<sys.npart;ipart++){
                        totbonds +=psl->pts[ipart].nbonds;
                    }
                    printf("tot bonds results from sum over all psl->pts[ipart].nbonds= %d \n", totbonds );
                    printf("nbonds before : %d  after : %d \n",nbonds_before, psl->nbonds);
                    error("ERROR:: nclusters not constant during clustermoves");
                }

                /*checks if any particle exeeds the maximum number of bonds (more than one per site)*/
                // printf("after cluster moves\n");
                check_maxbonds(psl);
            }
            else{
                total_energy(psl);
                // printf(" single moves, starting iwth %d bonds \n ",psl->nbonds);

                propagate_mc(psl);   
                // printf("a single particle move is performed, the energy is %lf\n", psl->energy);
                if (energy_divergence_check(psl,"single particle moves")==1) error("stop at single particle moves");
                /*checks if any particle exeeds the maximum number of bonds (more than one per site)*/
                
                check_maxbonds(psl); 
            }
        }
    }
    else{
        for(istep=0;istep<10;istep++){

            propagate_mc(psl);   

            if (energy_divergence_check(psl,"single particle moves")==1)   error("stop at single particle moves");

            /*checks if any particle exeeds the maximum number of bonds (more than one per site)*/
            check_maxbonds(psl);       
        } 
         
    }

    /*checks if any particle exeeds the maximum number of bonds (more than one per site)*/
    // printf("after single particle moves\n");
    check_maxbonds(psl);
    conf_output(psl);
     
    return ;
}

void propagate_mc(Slice *psl) {
    int i,nnlupdate=0, ipart;
    double movetype,dr2;
    int k,j;
    double r;


    //printf("propagating via mc\n");
    for(i=0; i<psl->nparts; i++) {
        movetype =  RandomIntegerRange(0, 2);
        single_particle_move(psl, movetype);
    }
    // energy_divergence_check(psl, "propagate_mc singlemoves");

    if ((analysis.xy_print==1) && (analysis.bond_breakage==0)) {
        // if you want to measure the bending ridigity of a chain
        // NOTE:    code doesn't check if the input is a chain. 
        //          this an easily be done setting start_type=2 in the input
        tailflipping(psl);
        
        // energy_divergence_check(psl, "propagate_mc tailflipping");
    }

    return;
}

void cluster_propagate_mc(Slice *psl){
    /* do a cluster move instead of single particle move*/
    int n, icluster;
    double nwhich;

    //printf("propagating via mc\n");
    for(n=0; n<psl->nparts; n++) {
        nwhich = RandomNumber();
        if(nwhich<.5) {
            icluster = translatepart_cluster_Ekparts(psl);
        }
        else{
            icluster = rotatepart_cluster_Ekparts(psl);
        }
    }

    return;
}


/* _______________ MOVEING SINGLE PARTICLES __________________________*/
void single_particle_move(Slice *psl, int r_or_t){
    /*one function for either rotate (1) or translate(0), selected via r_or_t*/
    int ipart, nbonds_old_tot,n,particle;
    Pts oldpart;
    double dr2, Eold, Enew, Edif;
    double tracked_diff, Eold_tot, Enew_tot;
    vector dr=nulvec;
    quaternion dq,dqrot;

    /*choose a random particle*/
    ipart=RandomIntegerRange(0, psl->nparts); 
    
    // old single particle info; and save the old particle info
    Eold=particle_energy(psl,ipart,0); 
    memcpy(&oldpart, &psl->pts[ipart], sizeof(Pts));

    if(r_or_t==0){
        // printf("translation \n ");
        mc_single.trans.tries++;
        /*make a random dr vector*/
        dr=RandomVector(mc_single.drmax);
        if (sys.gravity>0){
            dr.z/=(sys.gravity*10.);
        }
        /*displace the particle*/
        vector_add(psl->pts[ipart].r, dr, psl->pts[ipart].r);

        /*check if the particle is put into the wall, than always reject. No need to calculate the energy.*/
        if((sys.gravity>0) && (oldpart.r.z>1.) ){
            if(particle_in_wall(psl,ipart)==1){
                memcpy( &psl->pts[ipart],&oldpart, sizeof(Pts));

                return;
            }
        }
        // pbc(psl->pts[ipart].r,sys.boxl);    
    }
    else{
        // printf("rotation  ");
        mc_single.rot.tries++;

        /*rotate the quaternion*/
        dq=RandomQuaternionRange(mc_single.dqmax);
        rotate_quaternion_ipart(  psl,  ipart,  dq );
    }

    /*calculate the new particle energy*/
    Enew=particle_energy(psl,ipart,0);
    Edif = Enew - Eold;


    if (analysis.bond_breakage!=1){
        // reject too if bond breakage is not allowed
        if ( oldpart.nbonds != psl->pts[ipart].nbonds){
            /* if the translation caused a bond to break/form,if so do a clusterupdate before performing clustermoves*/
            memcpy( &psl->pts[ipart],&oldpart, sizeof(Pts));
            return;
        }
        else{
            for(n=0;n<oldpart.nbonds; n++){
                if(oldpart.bonds[n]!=psl->pts[ipart].bonds[n]){
                    memcpy( &psl->pts[ipart],&oldpart, sizeof(Pts));
                    return;
                }
            }
        }
    }

    // printf("%d  %.5lf  %.5lf  %.5lf \n",ipart, Eold , Edif,Enew);
    /*MC metropolis rule*/
    if(exp(-psl->beta*Edif)<RandomNumber()) {
        memcpy( &psl->pts[ipart],&oldpart, sizeof(Pts));
        return;
    }

    /*accepted*/
    psl->energy+=Edif;

    /*add acceptance count*/
    if (r_or_t==0){     mc_single.trans.acc++;}
    else{               mc_single.rot.acc++;}

    /* if the translation/rotation caused a bond to break/form, if so do a clusterupdate before performing clustermoves*/
    if (sys.cluster_MC ){
        if ( oldpart.nbonds != psl->pts[ipart].nbonds){ 
            cluster.update=1;
        }
        else{
            for(n=0;n<oldpart.nbonds; n++){
                if(oldpart.bonds[n]!=psl->pts[ipart].bonds[n]){
                    cluster.update=1;
                }
            }
        }
    }
    
    return;
}

/*_________________TAIL FLIPPING OF CHAIN_____________________________*/

void tailflipping(Slice *psl){
    /* choose a random colloid in the chain and flip the tail around the axis */
    quaternion dq, dqrot;
    tensor R;
    vector dr, rotvec;
    
    int ipart_ref, nclusters=psl->nclusters, icluster, ref_random;
    int  nbonds_old=0, nbonds_new=0, tot_nbonds;
    int ipart,jpart,i,n,direction;

    double Eold=0., Enew=0.,  Edif;
    double total_energy_old=psl->energy;
    Pts p_ref;

    /*select a cluster randomly*/
    icluster=(int)(RandomNumber()*nclusters);
    int clusteri_size=cluster.clustersizes[icluster];
    int cluster_nbonds = clusteri_size-1;
    Picl pici=cluster.pic[icluster];

    // /* only multi particle clusters */
    if(clusteri_size==1){
        return;
    }

    /*save old energy, nbonds, positions*/ 
    for(i=0; i<clusteri_size; i++){
        ipart=pici.stack[i];
        Eold+=particle_energy(psl,  ipart, 0)  ;
        nbonds_old+=(psl->pts[ipart].nbonds);
    }

    // the copy is used to make the new coordinates
    memcpy(cp_slice, psl, sizeof(Slice));

    /*pick randomly a reference particle*/
    ref_random= (int)(RandomNumber()*clusteri_size);
    ipart_ref = pici.stack[ref_random];
    p_ref = psl->pts[ipart_ref];


    /*the rotation quaternion and matrix based on the patchvector of the ref particle*/
    dq=RotateQuaternion(p_ref.patchvector[0],PI); 
    R = getrotmatrix(dq); 
    mc_tailflip.rot.tries++; 

    /*rotate the tail of the chain*/
    // step 1, make a list of all the particles in the chain 
    IntArray particlesleft_list; //

    //initiate particlesleft_list (filled)
    initIntArray(&particlesleft_list,  (size_t) clusteri_size);
    for(i=0; i<clusteri_size; i++){         
        insertIntArray(&particlesleft_list, pici.stack[i]);
    }
    // remove the reference particle
    removeElementXIntArray( &particlesleft_list ,  ipart_ref);

    // step 2, walk over the bonds till you reach the end
    // walk untill you reach particle with 1 bond or run into a particles thats not in the list. 
    direction=RandomIntegerRange(0, psl->pts[ipart_ref].nbonds); // randomly select forward or backward
    jpart=psl->pts[ipart_ref].bonds[direction]; // the first next particle


    if (checkElementXIntArray(&particlesleft_list, jpart )==1){
        // perform the tailflip here: 
        dr=particles_vector(psl,  ipart_ref, jpart); //return vector from ipart to jpart
        matrix_x_vector(R, dr, rotvec);
        vector_add(rotvec, cp_slice->pts[ipart_ref].r, cp_slice->pts[jpart].r)
    
        /* rotate the quaternion*/
        quat_times(dq,cp_slice->pts[jpart].q,dqrot);
        cp_slice->pts[jpart].q = dqrot;
        update_patch_vector_ipart(cp_slice,jpart); 

        DFS_tailflip( psl,jpart,  &particlesleft_list ,cp_slice , R, dq);
    }

    if(sys.gravity>0){
        for(i=0; i<cluster.clustersizes[icluster]; i++){
            ipart= pici.stack[i];
            if((psl->pts[ipart].r.z>=1.) && (particle_in_wall(cp_slice,ipart)==1)){
                freeIntArray(&particlesleft_list);
                return ;   
            }
        }
    }

    /*calc new  energy, nbonds*/
    /* perform energy caluclation always without using the neighborlist*/
    for(i=0; i<clusteri_size; i++){
        ipart = pici.stack[i];
        Enew+=particle_energy(cp_slice,  ipart, 0)  ;
        nbonds_new+=(cp_slice->pts[ipart].nbonds);
    }

    // potential energy coming from the bonds is not allowed to change during the 
    if (check_internal_energy( psl, cp_slice,  icluster, "tail flip")){
        error(" ERROR tail flip gone wrong. in bond energy\n");
        return;
    }

    Edif = (Enew - Eold);

    /*reject if*/
    if(exp(-sys.beta*Edif)<RandomNumber() || nbonds_old!=nbonds_new) {
        // printf(".     **reject tailflip Edif=%lf **\n",Edif);
        freeIntArray(&particlesleft_list);
        return ;
    }

    // printf(".     **accepted tailflip Edif=%lf  Enew=%lf Eold=%lf **\n",Edif, Enew, Eold);

    // free the memory
    freeIntArray(&particlesleft_list);
    // if accepted, copy the new positions to psl
    memcpy(psl,cp_slice, sizeof(Slice));   

    /* acccepted*/
    psl->energy=total_energy_old+Edif; 
    mc_tailflip.rot.acc++;
 
    return ;

}

void DFS_tailflip(Slice *psl,int ipart, IntArray *particlesleft_list, Slice *cp_slice , tensor R, quaternion dq){
    //  Depth First Search for tailflip
    vector dr;
    int jpart,j_inlist;
    vector rotvec;
    quaternion dqrot;

    // remove ipart from list, else you will find the bond back to ipart. Now you will walk forward to next bond
    removeElementXIntArray(particlesleft_list,  ipart );
    // dprint(ipart);
    // loop over the bound particles of ipart, to see if you have visited them
    for (int n = 0; n < psl->pts[ipart].nbonds; n++){   
        jpart=psl->pts[ipart].bonds[n]; // the bound particles to ipart;
  
        //  did you already visit jpart?
        j_inlist=checkElementXIntArray(particlesleft_list, jpart );
        if(j_inlist==1){
            
            // perform the tailflip here: 
            dr=particles_vector(psl,  ipart, jpart); //return vector from ipart to jpart
            matrix_x_vector(R, dr, rotvec);

            vector_add(rotvec, cp_slice->pts[ipart].r, cp_slice->pts[jpart].r)
   

            /* rotate the quaternion*/
            quat_times(dq,cp_slice->pts[jpart].q,dqrot);
            cp_slice->pts[jpart].q = dqrot;
            update_patch_vector_ipart(cp_slice,jpart); 

            DFS_tailflip( psl, jpart,  particlesleft_list,  cp_slice , R, dq );
        }
    }

    return;
}
/* _______________ MOVEING CLUSTERS __________________________________*/


int translatepart_cluster_Ekparts(Slice *psl) {

    int icluster, nclusters=psl->nclusters;
    int ipart, i, n;
    int update=0, nbonds_old=0, nbonds_new=0, tot_nbonds;
    double Ebond_new,Ebond_old;
    double Eold=0., Enew=0., Edif;
    vector dr;

    /*select a cluster randomly*/
    // icluster=(int)(RandomNumber()*nclusters);
    icluster = RandomIntegerRange(0,nclusters);
    if (icluster==nclusters) error("icluster==nclusters in translate cluster ");
    int clusteri_size=cluster.clustersizes[icluster];

    /* differentiate between single and multi particle clusters */
    if(clusteri_size==1){
        dr=RandomVector(mc_mono.drmax);
        mc_mono.trans.tries++;
    }
    else{ // printf("picked a cluster\n");
        dr=RandomVector(mc_cluster.drmax);
        mc_cluster.trans.tries++; 
    }
    if (sys.gravity>0){
        dr.z=0;
    }

    /*save old energy, nbonds, positions*/
    for(i=0; i<clusteri_size; i++){
        ipart=cluster.pic[icluster].stack[i];
        Eold+=particle_energy(psl, ipart,0) ;
        nbonds_old+=(psl->pts[ipart].nbonds);
    }
    
    /* save oldparticles; (!) might want to make this optimized such that you don't copy the whole structure*/
    memcpy(cp_slice , psl, sizeof(Slice));

    /*perform translation*/
     for(i=0; i<clusteri_size; i++){
        ipart = cluster.pic[icluster].stack[i];
        vector_add(psl->pts[ipart].r,dr,psl->pts[ipart].r); 
    } 

    /*check if any particle is put into the wall (z<1.0), than always reject. No need to calculate the energy.*/
    if(sys.gravity>0){
        for(i=0; i<cluster.clustersizes[icluster]; i++){
            ipart= cluster.pic[icluster].stack[i];
            if(cp_slice->pts[ipart].r.z>=1.){
                if(particle_in_wall(psl,ipart)==1){
                    memcpy(psl, cp_slice, sizeof(Slice));
                    // copy_clusterparticles( psl, cp_slice,icluster);
                    return -1;
                }
            }
        }
    }

    /*calc new  energy, nbonds*/
    /* perform energy caluclation always without using the neighborlist*/
    for(i=0; i<clusteri_size; i++){
        ipart = cluster.pic[icluster].stack[i];
        Enew+=particle_energy(psl, ipart,0) ;
        nbonds_new+=(psl->pts[ipart].nbonds);
    }

    //first check if you created bonds. then already reject due to detailed balance
    if (nbonds_new>nbonds_old){
        memcpy(psl, cp_slice, sizeof(Slice));
        // printf("rejected based on  nbonds\n");
        return -1;
    }
    else if (nbonds_new<nbonds_old){
        dprint(nbonds_new);
        dprint(nbonds_old);
        error("a bond has been broken during trnalstaional cluster move");
    }
    

     /* check if internal energy has stayed equal*/
    int trans_error= check_internal_energy( psl, cp_slice,  icluster, "cluster translation");

    if (trans_error){
        error("cluster translation gone wrong. in bond energy\n");
        // copy_clusterparticles( psl, cp_slice,icluster);
        return -1;
    }

    Edif = (Enew - Eold);

    /*if reject based on energy or nbonds */
    if((exp(-psl->beta*Edif)<RandomNumber() )|| (nbonds_old!=nbonds_new)) {
        memcpy(psl, cp_slice, sizeof(Slice));
        return -1;
    }

    psl->energy+=Edif;

    if(clusteri_size==1){   mc_mono.trans.acc++; }
    else{                   mc_cluster.trans.acc++;  }

    return icluster;
}

int rotatepart_cluster_Ekparts(Slice *psl) {

    int ipart_ref, nclusters=psl->nclusters, icluster, ref_random;
    int  nbonds_old=0, nbonds_new=0, tot_nbonds;
    int ipart,i,n;
    double Eold=0., Enew=0., Edif;
    double Ebond_new, Ebond_old;

    quaternion dq, dqrot;
    tensor R;
    vector dr, rotvec;

    /*select a cluster randomly*/
    icluster=(int)(RandomNumber()*nclusters);
    if (icluster==nclusters) error("icluster==nclusters in rotate cluster ");
    int clusteri_size=cluster.clustersizes[icluster];


    if(clusteri_size==1){
        rotate_monocluster(psl,icluster);
        return icluster;
    }
    else{ 
        if (sys.gravity>0){
            // if gravity , rotate around Z axis to prevent particles from being put inside wall
            dq=QuaternionZaxis((2*RandomNumber()-1)*mc_cluster.dqmax);  
        }
        else{
            dq=RandomQuaternionRange(mc_cluster.dqmax);
        }
        R = getrotmatrix(dq); 
        mc_cluster.rot.tries++; 
    }

    /*save old energy, nbonds*/
    for(i=0; i<clusteri_size; i++){
        ipart      = cluster.pic[icluster].stack[i];
        Eold       +=particle_energy(psl, ipart,0) ;
        nbonds_old +=(psl->pts[ipart].nbonds);
    }

    /*save old particles information (bv position, bond energy etc.)*/
    // copy_clusterparticles(  cp_slice,psl,icluster);
    memcpy(cp_slice,psl, sizeof(Slice));

    /*pick randomly the reference particle" ref_random -> [0,0.999> pick particle 0 etc*/
    ref_random = (int)(RandomNumber()*clusteri_size);
    Picl pici=cluster.pic[icluster]; 
    ipart_ref  = pici.stack[ref_random];
    Pts p_ref  = psl->pts[ipart_ref];
    vector dummy,dummy2;

    // printf("now rotating the structure\n");
    /*rotate the cluster, now you can just use the dr as defined below*/
    for(i=0; i<clusteri_size; i++){
        ipart= pici.stack[i];

        // dummy points from ipart_ref to ipart
        vector_minus(psl->pts[ipart].r,psl->pts[ipart_ref].r,dummy); 
        // no need to do pbc, because in linked_structure_icluster you already adjusted the coordinates 

        /* rotate the structure around the reference particle*/
        matrix_x_vector(R, dummy, dummy2);
        vector_add(dummy2,psl->pts[ipart_ref].r,psl->pts[ipart].r); // dummy points from ipart_ref to ipart

        /* rotate the quaternion and update patch vectors*/
        rotate_quaternion_ipart(  psl,  ipart,  dq );
    }

    /*calc new  energy, nbonds*/
    /* perform energy caluclation always without using the neighborlist*/
    for(i=0; i<clusteri_size; i++){
        ipart = pici.stack[i];
        Enew+=particle_energy(psl, ipart,0) ;
        nbonds_new+=(psl->pts[ipart].nbonds);
    }

    //first check if you created bonds. then already reject due to detailed balance
    if (nbonds_new>nbonds_old){
        memcpy(psl, cp_slice, sizeof(Slice));
        // copy_clusterparticles( psl, cp_slice,icluster);
        // printf("rejected based on  nbonds\n");
        return -1;
    }
    else if (nbonds_new<nbonds_old){
        error("WARNING: a bond has been broken during rotation cluster move");
    }

    /* check if internal energy has stayed equal*/
    if (check_internal_energy( psl, cp_slice,  icluster, " cluster rotation")){
        error(" WARNING cluster rotation gone wrong. in bond energy\n");
        return -1;
    }
    
    Edif = Enew - Eold;

    /*reject if*/
    if(exp(-psl->beta*Edif)<RandomNumber() ) {
        memcpy(psl, cp_slice, sizeof(Slice));
        // printf("rejected based on MC or nbonds\n");
        return -1;
    }

    /* acccepted*/
    psl->energy+=Edif;
    
    // printf("*** ROTATION SUCCESFULL **\n");
    if(clusteri_size!=1){
        mc_cluster.rot.acc++; 
    }
    
    return icluster;
}


/* _______________ MOVEING MONO CLUSTERS ______________________________*/
int rotate_monocluster(Slice *psl, int icluster){        

    quaternion dq, dqrot;
    tensor R;
    double Eold,Enew,Edif;
    int nbonds_new,nbonds_old;

    /*the rotation quaternion and matrix*/
    dq=RandomQuaternionRange(mc_mono.dqmax); 
    R = getrotmatrix(dq); 
    mc_mono.rot.tries++;
    int ipart= cluster.pic[icluster].stack[0];

    //old energy
    Eold=particle_energy(psl, ipart,0) ;
    nbonds_old=(psl->pts[ipart].nbonds);

    Pts oldpart =psl->pts[ipart];

    //the rotation
    rotate_quaternion_ipart(  psl,  ipart,  dq );

    //new energy + nbonds
    Enew=particle_energy(psl, ipart,0) ;
    nbonds_new=psl->pts[ipart].nbonds;

    Edif = (Enew - Eold);

    /*reject if*/
    if(exp(-psl->beta*Edif)<RandomNumber() || nbonds_old<nbonds_new) {
        // memcpy(psl->pts, oldpart, sizeof(psl->pts));
        psl->pts[ipart]=oldpart;
        // printf("rejected based on MC or nbonds\n");
        return -1;
    }

    /* acccepted*/
    psl->energy+=Edif;
    
    mc_mono.rot.acc++;

    return icluster;
}

void copy_clusterparticles( Slice *psl_new ,Slice *psl_source,int icluster){
    // copies the particles in icluster from psl_source to psl_new
    int clusteri_size=cluster.clustersizes[icluster];
    int ipart;

    /*perform copy on selected particles */
    for(int i=0; i<clusteri_size; i++){
        ipart = cluster.pic[icluster].stack[i];
        psl_new->pts[ipart]=psl_source->pts[ipart];
    } 


    return;
}



/* _______________ CHECK, OPTIMIZE AND PRINT MC PARAMETERS _____________*/
int energy_divergence_check(Slice *psl,  char loc[50]){
    // checks for the divergence of the energy due to the use of a running energy in the MC code
    double diff = fabs(total_energy(psl)-psl->energy);
    // printf("chek energy_divergence_check\n");
    if(diff>.001){
        printf("\nWARINING: there is a difference of %.6lf in the recalc.  after %s\n", diff,loc);
        printf("recalc total_energy              %lf\n", total_energy(psl));
        printf("tracked slice[0].energy               %lf\n", psl->energy);
        print_slice_information(psl);
        if(diff>1){
            printf("the energy difference is too big\n ");
            // print_slice_information(psl);
            error("stop");
        }
        return 1;
    }
    return -1;
}

int check_internal_energy(Slice *psl_new,Slice *psl_old, int icluster, char movetype[100]){
    int clusteri_size=cluster.clustersizes[icluster];
    double Ebond_new,Ebond_old;
    int cluster_error=0;

    for(int i=0; i<clusteri_size; i++){   // loop over cluster 
        int ipart= cluster.pic[icluster].stack[i];
        for(int n=0;n<psl_new->pts[ipart].nbonds;n++){ // loop over particle in the cluster and look at bonds

            // 
            if ( psl_old->pts[ipart].nbonds != psl_new->pts[ipart].nbonds){ // "<"bonds are broken? 

                printf(" in %s  number of bonds has changed from %d (old) to %d (new) for particle %d\n", movetype,psl_old->pts[ipart].nbonds,psl_new->pts[ipart].nbonds, ipart);
                // continue;
                for(int m=0;m<psl_new->pts[ipart].nbonds;m++){
            
                    Ebond_new=psl_new->pts[ipart].bond_energy[m];
                    Ebond_old=psl_old->pts[ipart].bond_energy[m];
                    dprint(m);
                    dprint(psl_new->pts[ipart].bonds[m]);
                    dprint(psl_old->pts[ipart].bonds[m]);
                    if(fabs(Ebond_new-Ebond_old) >1e-9){
                        printf(" ipart: %d,  Ebond_new = %.12lf , Ebond_old=%.12lf,      difference = %.12f\n",ipart,Ebond_new,Ebond_old, fabs(Ebond_new-Ebond_old));
                        
                        cluster_error=1;

                    }

                    printf("    psl the other particle     %d has %d bonds \n",psl_new->pts[ipart].bonds[m],psl_new->pts[psl_new->pts[ipart].bonds[m]].nbonds);
                    printf("    psl_old the other particle %d has %d bonds \n",psl_old->pts[ipart].bonds[m],psl_old->pts[psl_old->pts[ipart].bonds[m]].nbonds);
                }
                printf("\n");
                for(int m=0;m<psl_old->pts[ipart].nbonds;m++){
            
                    Ebond_new=psl_new->pts[ipart].bond_energy[m];
                    Ebond_old=psl_old->pts[ipart].bond_energy[m];
                    dprint(m);
                    dprint(psl_new->pts[ipart].bonds[m]);
                    dprint(psl_old->pts[ipart].bonds[m]);
                    if(fabs(Ebond_new-Ebond_old) >1e-9){
                                        
                        
                        printf("ipart: %d, Ebond_new = %.12lf , Ebond_old=%.12lf,      difference = %.12f\n",ipart,Ebond_new,Ebond_old, fabs(Ebond_new-Ebond_old));
                        
                        cluster_error=1;
                        // printf("rejected based on internal bond_energy\n");
                        // return -1;
                    }

                    printf(" psl the other particle %d has %d bonds \n",psl_new->pts[ipart].bonds[m],psl_new->pts[psl_new->pts[ipart].bonds[m]].nbonds);
                    printf(" psl_old the other particle %d has %d bonds \n",psl_old->pts[ipart].bonds[m],psl_old->pts[psl_old->pts[ipart].bonds[m]].nbonds);
                }
                printf("\n");
            }
        }
    }
    return cluster_error;
}




void printstatusmc_sub(MC  *type){
        
    if ((type->rot.tries>0) || (type->trans.tries>0))  printf("  %s acceptance:\n",type->name);
    
    /*calculate the ratio*/
    type->rot.ratio=(double)type->rot.acc/(double)type->rot.tries;
    type->trans.ratio=(double)type->trans.acc/(double)type->trans.tries;

    if ((type->trans.tries>0))printf("   trans:    %6ld/%6ld   = %.3lf   >  drmax = %.3lf\n",type->trans.acc,type->trans.tries,type->trans.ratio,type->drmax);
    if ((type->rot.tries>0) )printf("   rot:      %6ld/%6ld   = %.3lf   >  dqmax = %.3lf\n",type->rot.acc,type->rot.tries,type->rot.ratio,type->dqmax);

    return;
}


void printstatusmc() {
    static int init=1;
    double en,l;
    
    if ( init ){
        init=0;
        return;
    }
    
    printstatusmc_sub(&mc_single);
    if(sys.cluster_MC==1){
        printstatusmc_sub(&mc_cluster);
        printstatusmc_sub(&mc_mono);
    }

    if ((analysis.xy_print==1) && (analysis.bond_breakage==0)){
        printstatusmc_sub(&mc_tailflip);
    }

    return;
}

void optimizemc_sub(MC *mctype) {

    static int initiate=1;

    if (initiate ){
        initiate=0;
        return;
    }
    /*calculate the ratio*/
    mctype->rot.ratio=(double)mctype->rot.acc/(double)mctype->rot.tries;
    mctype->trans.ratio=(double)mctype->trans.acc/(double)mctype->trans.tries;

    /*update the total accepted, tried translate or rotate moves. */
    mctype->finrot.acc+=mctype->rot.acc;
    mctype->finrot.tries+=mctype->rot.tries;
    mctype->fintrans.acc+=mctype->trans.acc;
    mctype->fintrans.tries+=mctype->trans.tries;
    
    
   
    // update the max translation and rotation
    /*rotation*/
    if((mctype->rot.ratio<0.3 )&& (mctype->dqmax>0.0001)) {
        mctype->dqmax/=1.1;
    }
    else if((mctype->rot.ratio>0.7) && (mctype->dqmax<179)) {
        mctype->dqmax*=1.1;
    }
    else if(mctype->dqmax>180){
        mctype->dqmax=180.;
    }

    /*translation*/
    if((mctype->trans.ratio<0.3) && (mctype->drmax>0.0001)) {
        mctype->drmax/=1.1;
    }
    else if((mctype->trans.ratio>0.7) && (mctype->drmax<=sys.boxl.x)) {
        mctype->drmax*=1.1;
    }
    else if (mctype->drmax>sys.boxl.x){
        mctype->drmax=sys.boxl.x;
    }

    /*set accept. and tries to zero*/
    mctype->rot.acc=0;
    mctype->rot.tries=0;

    mctype->trans.acc=0;
    mctype->trans.tries=0;
    

    return;
}

void optimizemc() {
    static int init=1;
    double en,l;
    
    if ( init ){
        init=0;
        return;
    }
   
    optimizemc_sub(&mc_single);

    if(sys.cluster_MC==1){
        optimizemc_sub(&mc_cluster);
        optimizemc_sub(&mc_mono);
    }

    return;
}



void final_printstatusmc_sub(MC type){
        
    printf("%s\n",type.name);
    type.fintrans.ratio=(double)type.fintrans.acc/(double)type.fintrans.tries;
    type.finrot.ratio=(double)type.finrot.acc/(double)type.finrot.tries;

    printf("Accepted %6ld translations from %6ld tries, ratio %lf\n",type.fintrans.acc,type.fintrans.tries,type.fintrans.ratio);
    printf("Accepted %6ld rotations    from %6ld tries, ratio %lf\n",type.finrot.acc,type.finrot.tries,type.finrot.ratio);

    return;
}

/* _______________ SETTING UP MC ________________________________________*/
void setup_MC(void){
    
    //Monte Carlo
    printf("\nSetting up the MC parameters... ");
    //Monte Carlo; there is only MC in this code...

    if ((sys.nearest_neighbor==1) & (sys.sim_type==MC_ALGORITHM) ){
        error("nearest_neighbor and MC cannot go together! turn nearest_neighbor off ");
    }

    cp_slice=(Slice *)calloc(1,sizeof(Slice)); //global variable,used in cluster MC
    copyslice=(Slice *)calloc(1,sizeof(Slice)); //global variable,used in cluster MC
    printf("\n **Setting up the MC parameters...\n ");

    
    sprintf(mc_single.name,"single particle");
    setup_mc_move(&mc_single);
    // mc_single.drmax/=10.;
    // mc_single.dqmax/=10.;

    
    if(sys.cluster_MC==1){
        sprintf(mc_cluster.name,"cluster moves");
        setup_mc_move(&mc_cluster);

        sprintf(mc_mono.name,"single-particle-clusters moves");
        setup_mc_move(&mc_mono);
    }  

    if ((analysis.xy_print==1) && (analysis.bond_breakage==0)){
        setup_mc_move(&mc_tailflip);
    }
    printf("   .. done **\n");  
    
    return;
}

void setup_mc_move(MC *mctype){
        printf("*****setting up the MC %s ****\n",mctype->name);
        mctype->rot.acc=0;
        mctype->trans.acc=0;
        mctype->fintrans.acc=0;
        mctype->finrot.acc=0;
        mctype->trans.tries=0;
        mctype->fintrans.tries=0;
        mctype->rot.tries=0;
        mctype->finrot.tries=0;

        mctype->drmax=1.;
        mctype->dqmax=5.;
        return;
}