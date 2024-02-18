#include "path.h"
/* _______________ THE MONTE CARLO ALGORITHM __________________________*/

void mccycle(Slice *psl) {
    /*sys.sim_type==2
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
            if(which<1.){
                total_energy(psl); // it's necessary to do total_energy before knowing nbonds 
                // printf(" clustermoves, starting iwth %d bonds \n ",psl->nbonds);
                nbonds_before= psl->nbonds;
                // copy cluster if things went wrong to go back in time.
                memcpy(copyslice,psl,sizeof(Slice));

                /* expensive function, do this function only when the number of bonds has changed. */
                if (cluster.update==1){
                    printf("****performing cluster analysis and linking*****\n");
                    linking_all_cluster(psl);
                    cluster.update=0;

                    printf("  there are %d clusters \n",psl->nclusters );
                }

                
                /* *** the cluster move*** */
                cluster_propagate_mc(psl);

                // printf("a cluster move is performed, the energy is         %lf\n", psl->energy);
                Edivcheck=energy_divergence_check(psl,"cluster move");
                if (Edivcheck>-1){
                    printf("something went wrong in cluster moves, to back to snapshot before cluster moves\n");

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
            total_energy(psl);
            // printf(" single moves, starting iwth %d bonds \n ",psl->nbonds);
            linking_all_cluster(psl);
            propagate_mc(psl);   
            // printf("a single particle move is performed, the energy is %lf\n", psl->energy);

            if (energy_divergence_check(psl,"single particle moves")==1) error("stop at single particle moves");
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
        ipart = single_particle_move(psl, movetype);
    }
    return;
}


/* _______________ MOVEING SINGLE PARTICLES __________________________*/
int single_particle_move(Slice *psl, int r_or_t){
    /*one function for either rotate (0) or translate(1)*/
    int ipart, nbonds_old_tot,n,particle;
    Pts oldpart;
    double dr2, Eold, Enew, Edif;
    double tracked_diff, Eold_tot, Enew_tot;
    vector dr, dr_dummy;
    quaternion dq,dqrot;

    /*choose a random particle*/
    ipart=(int)(RandomNumber()*psl->nparts);
    
    // old single particle info
    // printf("Eold %d :\n",ipart);
    // Eold_tot=total_energy(psl);
    Eold=particle_energy(psl,ipart,0); 
    oldpart=psl->pts[ipart];

    if((psl->energy/psl->energy!=1) && (fabs(psl->energy)>1e-10)){
        error("energy nan in before  single particle move ");
    }

    MC *mc_single; 
    // DETAILED BALANCE NOT OBEYED WITH THIS STRUCTURE. OR DONT USE OR ADJUST ACCEPTANCE RULE
    // if the particles has no bonds, give is mc_single_large , the drmax and dq ax are larger
    // if its bound, these parameres are much smaller
    // if (psl->pts[ipart].nbonds==0) mc_single=&mc_single_large ;
    // else mc_single=&mc_single_small ;

    mc_single=&mc_single_small;

    if (r_or_t==0){
        // printf("translation  ");
        mc_single->trans.tries++;
        /*make a random dr vector*/
        dr=RandomVector(mc_single->drmax);
        if (sys.gravity>0){
            dr.z/=(sys.gravity*10.);
        }
        /*displace the particle*/
        vector_add(psl->pts[ipart].r, dr, psl->pts[ipart].r);

        /*check if the particle is put into the wall, than always reject. No need to calculate the energy.*/
        if(sys.gravity>0){
            if(oldpart.r.z>1.){
            // if(oldpart[ipart].r.z>1.){
                if(particle_in_wall(psl,ipart)==1){
                    psl->pts[ipart]= oldpart;
                    return -1;
                }
            }
        }
        // pbc(psl->pts[ipart].r,sys.boxl); /* perform pbc if no neighborlist is used*/   
    }
    else{
        // printf("rotation  ");
        mc_single->rot.tries++;

        /*rotate the quaternion*/
        dq=RandomQuaternionRange(mc_single->dqmax);
        // dq=QuaternionZaxis( 1.);
        rotate_quaternion_ipart(  psl,  ipart,  dq );
    }

    /*calculate the new particle energy*/
    // printf("Enew %d :\n",ipart);
    // Enew_tot=total_energy(psl);
    Enew=particle_energy(psl,ipart,0);
    Edif = Enew - Eold;


    if (analysis.bond_breakage!=1){
        // reject too if bond breakage is not allowed
        if ( oldpart.nbonds != psl->pts[ipart].nbonds){
            /* if the translation caused a bond to break/form,if so do a clusterupdate before performing clustermoves*/
            psl->pts[ipart]=oldpart;
            return -1;
        }
        else{
            for(n=0;n<oldpart.nbonds; n++){
                if(oldpart.bonds[n]!=psl->pts[ipart].bonds[n]){
                    psl->pts[ipart]=oldpart;
                    return -1;
                }
            }
        }
    }
    // printf("%d  %.5lf  %.5lf  %.5lf \n",ipart, Eold , Edif,Enew);
    /*metropolis rule*/
    if(exp(-psl->beta*Edif)<RandomNumber()) {
        // tracked_diff=Eold_tot-psl->energy;
        // if (fabs(tracked_diff)>0.0001){
        //     gprint(Edif);
        //     printf("WARNING difference of %.10lf from Etot_old =  %.5lf and psl->energy= %.5lf\n",tracked_diff,Eold_tot,psl->energy);
        // }


        psl->pts[ipart]= oldpart;
        return -1;
    }


    /*accepted*/
    psl->energy+=Edif;

    if((psl->energy/psl->energy!=1) && (fabs(psl->energy)>1e-10)){
        printf("Edif %lf\n", Edif);
        printf("Enew %lf -Eold %lf\n",Enew,Eold );
        error("energy nan single particle move");
    }
    
    // // tracked_diff=Enew_tot-psl->energy;

    // if (fabs(tracked_diff)>0.0001){
    //     printf("WARNING difference of %.10lf from Etot =  %.10lf and psl->energy= %.10lf\n",tracked_diff,Enew_tot,psl->energy);
    // }
    // printf("accepted %d. Edif =  %.5lf  - %.5lf   = %.5lf      Etot=%.12lf\n\n",ipart, Enew, Eold,Edif,  psl->energy);

    /*add acceptance*/
    if (r_or_t==0){
        mc_single->trans.acc++;
        particle=-1;
    }
    else{
        mc_single->rot.acc++;
        particle = ipart; //return the particle number for updating dr (nearestneighbor)
    }

    /* if the translation/rotation caused a bond to break/form,if so do a clusterupdate before performing clustermoves*/
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
    // energy_divergence_check(psl,"single particle moves");
    // printf("\n");

    return particle;
}

/* _______________ MOVEING CLUSTERS __________________________________*/
void cluster_propagate_mc(Slice *psl){
    /* do a cluster move instead of single particle move*/

    int n,update=0, i, icluster, ipart;
    double nwhich,dr2;

    //printf("propagating via mc\n");
    for(n=0; n<psl->nparts; n++) {
        nwhich = RandomNumber();
        if(nwhich<1.) {
            icluster = translatepart_cluster_Ekparts(psl);
        }
        else{
            icluster = rotatepart_cluster_Ekparts(psl);
        }
    }

    return;
}

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
    int clusteri_size=cluster.clustersize[icluster];



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
        dr.z/=(sys.gravity*50.);
        dr.z=0;
    }
   

    /*save old energy, nbonds, positions*/
    for(i=0; i<clusteri_size; i++){
        ipart=cluster.pic[icluster].stack[i];
        Eold+=particle_energy(psl, ipart,0) ;
        nbonds_old+=(psl->pts[ipart].nbonds);
    }
   
    
    /* save oldparticles; (!) might want to make this optimized such that you don't copy the whole structure*/
    memcpy(psl_old , psl, sizeof(Slice));
    // copy_clusterparticles(psl_old , psl, icluster);

    /*perform translation*/
     for(i=0; i<clusteri_size; i++){
        ipart = cluster.pic[icluster].stack[i];
        vector_add(psl->pts[ipart].r,dr,psl->pts[ipart].r); 
    } 

    /*check if any particle is put into the wall (z<1.0), than always reject. No need to calculate the energy.*/
    if(sys.gravity>0){
        for(i=0; i<cluster.clustersize[icluster]; i++){
            ipart= cluster.pic[icluster].stack[i];
            if(psl_old->pts[ipart].r.z>=1.){
                if(particle_in_wall(psl,ipart)==1){
                    memcpy(psl, psl_old, sizeof(Slice));
                    // copy_clusterparticles( psl, psl_old,icluster);
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
        memcpy(psl, psl_old, sizeof(Slice));
        // copy_clusterparticles( psl, psl_old,icluster);
        // printf("rejected based on  nbonds\n");
        return -1;
    }
    else if (nbonds_new<nbonds_old){
        dprint(nbonds_new);
        dprint(nbonds_old);
        error("a bond has been broken during trnalstaional cluster move");
    }
    

     /* check if internal energy has stayed equal*/
    int trans_error= check_internal_energy( psl, psl_old,  icluster, "cluster translation");

    if (trans_error){
        error("cluster translation gone wrong. in bond energy\n");
        // copy_clusterparticles( psl, psl_old,icluster);
        return -1;
    }
        

    Edif = (Enew - Eold);

    /*if reject based on energy or nbonds */
    if((exp(-psl->beta*Edif)<RandomNumber() )|| (nbonds_old!=nbonds_new)) {
        memcpy(psl, psl_old, sizeof(Slice));
        // copy_clusterparticles( psl, psl_old,icluster);
        return -1;
    }

    psl->energy+=Edif;


    if(clusteri_size==1){
        mc_mono.trans.acc++;
    }
    else{
        mc_cluster.trans.acc++; 
    }

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
    // icluster=(int)(RandomNumber()*nclusters);
    if (icluster==nclusters) error("icluster==nclusters in rotate cluster ");
    int clusteri_size=cluster.clustersize[icluster];


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
    // copy_clusterparticles(  psl_old,psl,icluster);
    memcpy(psl_old,psl, sizeof(Slice));

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
        memcpy(psl, psl_old, sizeof(Slice));
        // copy_clusterparticles( psl, psl_old,icluster);
        // printf("rejected based on  nbonds\n");
        return -1;
    }
    else if (nbonds_new<nbonds_old){
        error("WARNING: a bond has been broken during rotation cluster move");
    }

    /* check if internal energy has stayed equal*/
    if (check_internal_energy( psl, psl_old,  icluster, " cluster rotation")){
        error(" WARNING cluster rotation gone wrong. in bond energy\n");
        return -1;
    }
    
    Edif = Enew - Eold;

    /*reject if*/
    if(exp(-psl->beta*Edif)<RandomNumber() ) {
        memcpy(psl, psl_old, sizeof(Slice));

        // copy_clusterparticles( psl, psl_old,icluster);

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
    int clusteri_size=cluster.clustersize[icluster];
    int ipart;

    /*perform copy on selected particles */
    for(int i=0; i<clusteri_size; i++){
        ipart = cluster.pic[icluster].stack[i];
        psl_new->pts[ipart]=psl_source->pts[ipart];
    } 


    return;
}



/* _______________ CHECK, OPTIMIZE AND PRINT MC PARAMETERS _____________*/

int energy_divergence_check(Slice *psl, char loc[50]){
    // checks for the divergence of the energy due to the use of a running energy in the MC code
    double diff = fabs(total_energy(psl)-psl->energy);
    // printf("chek energy_divergence_check\n");
    if(diff>.001){
        printf("WARINING: there is a difference of %.6lf in the recalc.  after %s\n", diff,loc);
        printf("recalc total_energy              %lf\n", total_energy(psl));
        printf("tracked slice[0].energy               %lf\n", psl->energy);
        psl->energy=total_energy(psl);
        if(diff>1){
            printf("there is still a difference of %.6lf  in the recalc.  after %s\n", diff, loc);
            printf("recalc total_energy              %lf\n", total_energy(psl));
            printf("tracked slice[0].energy               %lf\n", psl->energy);
            return 1;

        }
        return 0;

    }
    return -1;
}

int check_internal_energy(Slice *psl_new,Slice *psl_old, int icluster, char movetype[100]){
    int clusteri_size=cluster.clustersize[icluster];
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
    
    printstatusmc_sub(&mc_single_large);
    printstatusmc_sub(&mc_single_small);
    // printf("The address of mc_single is %p \n",  &mc_single);
    if(sys.cluster_MC==1){
        printstatusmc_sub(&mc_cluster);
        // printf("The address of mc_cluster is %p \n",  &mc_cluster);
        printstatusmc_sub(&mc_mono);
        // printf("The address of mc_mono is %p \n",  &mc_mono);

    }

    return;
}

void optimizemc_sub(MC *mctype) {

    /*single particle moves*/
    static int initiate=1;

    if (initiate ){
        initiate=0;
        return;
    }
    /*calculate the ratio*/
    mctype->rot.ratio=(double)mctype->rot.acc/(double)mctype->rot.tries;
    mctype->trans.ratio=(double)mctype->trans.acc/(double)mctype->trans.tries;

    /*update the total accepted, tried translate or rotate moves*/
    mctype->finrot.acc+=mctype->rot.acc;
    mctype->finrot.tries+=mctype->rot.tries;
    mctype->fintrans.acc+=mctype->trans.acc;
    mctype->fintrans.tries+=mctype->trans.tries;
    // printf("ratio = %d / %.d\n",rot.acc,rot.tries);
    
    /*set accept. and tries to zero*/
    mctype->rot.acc=0;
    mctype->rot.tries=0;

    mctype->trans.acc=0;
    mctype->trans.tries=0;
   
    // printf("%s\n", mctype->name);
    /*rotation*/
    if((mctype->rot.ratio<0.3 )&& (mctype->dqmax>0.0001)) {
        // printf(" rot decrease\n");
        // gprint(mctype->rot.ratio);
        mctype->dqmax/=1.1;
        // gprint(mctype->dqmax);
    }
    else if((mctype->rot.ratio>0.7) && (mctype->dqmax<179)) {
        // printf(" rot increase\n");
        // gprint(mctype->rot.ratio);

        mctype->dqmax*=1.1;
         // gprint(mctype->dqmax);
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
        // gprint(mctype->drmax);
    }
    else if (mctype->drmax>sys.boxl.x){
        mctype->drmax=sys.boxl.x;
    }
    

    return;
}

void optimizemc() {
    static int init=1;
    double en,l;
    
    if ( init ){
        init=0;
        return;
    }
   
    optimizemc_sub(&mc_single_large);
    optimizemc_sub(&mc_single_small);

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

    if (sys.nearest_neighbor==1 & sys.sim_type==MC_ALGORITHM ){
        error("nearest_neighbor and MC do not go together!");
    }

    psl_old=(Slice *)calloc(1,sizeof(Slice)); //global variable,used in cluster MC
    copyslice=(Slice *)calloc(1,sizeof(Slice)); //global variable,used in cluster MC
    printf("\nSetting up the MC parameters...\n ");

    sprintf(mc_single_large.name,"single particle large moves");
    setup_mc_move(&mc_single_large);

    sprintf(mc_single_small.name,"single particle small moves");
    setup_mc_move(&mc_single_small);
    mc_single_small.drmax/=100.;
    mc_single_small.dqmax/=10.;

    gprint(mc_single_large.drmax);
    gprint(mc_single_large.dqmax);
    
    if(sys.cluster_MC==1){
        sprintf(mc_cluster.name,"cluster moves");
        setup_mc_move(&mc_cluster);

        sprintf(mc_mono.name,"single-particle-clusters moves");
        setup_mc_move(&mc_mono);
    }  
    printf("   .. done\n");  
    
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