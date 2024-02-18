#include "path.h"

/*------------------LOCAL FUNCTIONS------------------------------------------*/
double E_active_alignment(Slice *, int );
void particle_pair_energy(Slice *, int , int , double *, double *, double *);
double total_Epair_neighborlist(Slice *);
double wall_interaction_ipart(Slice *, int );
double total_Epair_woneighborlist(Slice *);
/*---------------------------------------------------------------------------*/


double total_energy(Slice *psl) { 
    /*  this function calculates the total energy (= gravitational energy + bond energies) of the system.
    A loop over all particles is done.
    */
    double Etot=0,Egrav_tot=0,Ewall_tot=0,Eactiveforce=0,Ee2e=0,Epair=0;
    Pts *psi;
  
    /* gravity*/
    for( int ipart=0; ipart<psl->nparts; ipart++ )  {
        psl->pts[ipart].nbonds =0; // set all nbonds to zero, will be updated in total_Epair
        psi = &psl->pts[ipart];
        if(sys.gravity>0){
            Egrav_tot+=gravitational_energy_ipart(psl, ipart);
            if (pot.wall_int!=0){                           Ewall_tot+=wall_interaction_ipart(psl,ipart);  }
            if (sys.particletype[psi->ptype].activity==1){  Eactiveforce+=E_active_alignment(psl,  ipart); }
        }  
    }


    // Vpair with or without neighbor list, loops differently over the particles
    if (sys.nearest_neighbor){
        // perform energy calculation with nearest neighbor list
        // printf("with nn \n");
        Epair=total_Epair_neighborlist(psl);
    }
    else{
        // printf("without nn \n");
        Epair=total_Epair_woneighborlist(psl);
    }

    Etot=Epair+Egrav_tot+Ewall_tot+Eactiveforce+Ee2e;
    // printf("The total energy is %.5lf\n",Etot );
    return Etot;
}

vector find_directing_patchvector(Slice *psl, int ipart,vector u1){
    /* find the patchevector on particle ipart that makes the smalles  angle with interparticle vector u1
    it does track here if this angle falls within the patch range/width*/
    int isite; 
    double tracker_i=-1,cositheta;
    vector p1_select=nulvec;

    for( isite=0; isite<sys.particletype[psl->pts[ipart].ptype].nsites; isite++ ) {
        vector p1=psl->pts[ipart].patchvector[isite];
        cositheta=vector_inp(u1,p1); 
        
        if(cositheta<tracker_i) {
            continue;
        }
        tracker_i=cositheta;
        p1_select=p1;
    }

    return p1_select;
}

double harmonic_oscillator(double k,  double x_shift,double x){
    return 1./2*k*(x-x_shift)*(x-x_shift);
}

double wall_interaction_ipart(Slice *psl, int ipart){
    /* if the wall has an interaction with the patches
    there is a harmonic oscillator that governs the interaction strength with its minimum below 0*/
    Pts *psi;
    psi=&psl->pts[ipart];
    int isite;
    vector z_axis=nulvec,x_axis=nulvec,y_axis=nulvec,p1;
    double Ewall_tmp=0,Ewall,z,angle_max=site[psi->ptype].cosdelta;
    double zcut=sys.particletype[psi->ptype].zcut+sys.particletype[psi->ptype].radius*1.1,angle,cosangle,S,E_height,pi180=PI/180.;
    z_axis.z=1;
    x_axis.x=1;
    y_axis.y=1;

   
    // only if the particle is in the vicinity of the wall: 
    if (psi->r.z<=zcut & psi->r.z>0){
        // find thenpatch that is pointing downwards. take the innerproduct with patch vector if 
        for(isite=0;isite<sys.particletype[psi->ptype].nsites;isite++){
            p1=psi->patchvector[isite];
            if (p1.z>0){
                continue;
            }
            cosangle=-vector_inp(p1,z_axis); 
            
            if (cosangle<angle_max){
                continue;
            }
            // printf("wall angle ipart %d = %lf\n", ipart,angle);
            if (psi->ptype>PTYPES){
                dprint(PTYPES);
                dprint(psi->ptype);
                error("the Ptype is too large");
            }
            S= Saccent( cosangle, psi->ptype);
            
            // E_height=gravitational_energy_ipart(psl, ipart); // gravitational energy is zero at the wall
            // z=fabs(p1.z-zcut);
            E_height=-1.*harmonic_oscillator(pot.wall_int,zcut,p1.z);  
            // E_height=0;
            Ewall=E_height*S;

            // gprint(Ewall);
            if (Ewall<Ewall_tmp){
                Ewall_tmp=Ewall;

                printf("*** ipart %d with isite %d***\n",ipart,isite);
                gprint(Ewall);
                gprint(S);
                angle=cosangle_to_angle(cosangle);
                gprint(angle);
                printf("\n");

            }
        }
    }
    
    return Ewall_tmp ;//-Ewall_max;
}


double E_active_alignment(Slice *psl, int ipart){
    // the energy penalty for misalignment with wall of active force
    Pts *psi;
    psi = &psl->pts[ipart];

    double cosangle,angle,rad_to_degree=180./PI,Eactiveforce;
    vector Fa;
    
    tensor rotmat = getrotmatrix(psi->q); 
    matrix_x_vector(rotmat,sys.particletype[psi->ptype].e_A,Fa); // e_A unit vector; Fa the direction of the active force
   
    angle=fabs(asin(Fa.z));

    Eactiveforce=0.5*500*angle*angle; // used to be 500=0.5*1000

    return Eactiveforce;
}

double particle_energy(Slice *psl, int ipart, int jpart) { 
    /*  this function calculates the particle energy of particle ipart.
    this includes gravity is on.
    It starts counting from particle jpart.
    Use jpart=0 in MC code, while in total_energy you use jpart=ipart+1;

    For MC you need to count bonds, 
    for the calculation of a single particle you need to known if nr bonds changed (fot the cluster update)
    Only info about that specific particle is needed

    If total energy is calculated, all the particles should contain info about bonden (so if p1-p2 bond, than both should know that)
    */

    Pts *psi ;    
    double Ebond=0,Erep=0;
    double Ebond_tot=0, Etot=0, Erep_tot=0, Egrav=0, Ewall=0,Eactiveforce=0,Ee2e=0;
    double s;
    int kpart;
    // printf("  particle_energy \n");
    
    psi = &psl->pts[ipart];

    if (jpart==0){
        //if you start from MC (jpart==0) then set nbonds of ipart to zero first
        //if you start from total_energy it is already set to zero in that function
        psl->pts[ipart].nbonds =0;

        // if end-to-end restrained
        //  && only do when coming from MC (if jpart==0), else (if jpart>0) , ie. via total_energy), it's done outside particle_energy().
    
    }   

    /* gravity*/
    if(sys.gravity>0){
        Egrav=gravitational_energy_ipart(psl, ipart);
        if (pot.wall_int!=0){  Ewall=wall_interaction_ipart(psl,ipart);     }
        if (sys.particletype[psi->ptype].activity==1){  Eactiveforce=E_active_alignment(psl,  ipart); }
    }

    /*kpart is de teller, bereken de energy van ipart, start counting the sum vanaf jpart*/
    for( kpart=jpart; kpart<sys.npart; kpart++ ) {
  
        // skip jezelf
        if(ipart==kpart){   continue; }


        particle_pair_energy(psl,ipart,kpart,&Erep, &Ebond, &s);


        Erep_tot += Erep; 
        Ebond_tot += Ebond ;
        if(Ebond<pot.bond_cutoffE){
            add_bond_information(psl,ipart,kpart,Erep,Ebond,s);
            
            /*in case of calculating the total energy, also add the bond to the other particle*/
            if( (jpart>0 )&& (kpart>ipart)){      
                add_bond_information(psl,kpart,ipart,Erep,Ebond,s);
            }        
        }
    }

    Etot = Ebond_tot + Erep_tot + Egrav + Ewall + Eactiveforce + Ee2e;
    // printf("The energy of particle %d is %.5lf =  %.5lf (Ebond_tot) +  %.5lf (Erep_tot) +  %.5lf (Egrav)  +  %.5lf (Ewall) +  %.5lf (Eactiveforce)+  %.5lf (Ee2e)\n",ipart,Etot,Ebond_tot,Erep_tot ,Egrav ,Ewall, Eactiveforce , Ee2e);
    return Etot;
}

void orientation_parameters(Slice *psl, int ipart, int jpart, vector u1, double *costheta_i, double *costheta_j, double *costheta_ij){
    /* this function calculates the bulk0bulk distance s, the bulk-patch distance s_bp, the patch-pathc distance s_pp,
    the theta and phi angles */

    Pts *psi,*psj;
    
    vector min_u1, p1,p2;
    int calc_ij=0; // this is 0 if you dont use S180, hardcoded?

    psi=&psl->pts[ipart];
    psj=&psl->pts[jpart];

    //find smallest combination of the patch vector angles
    scalar_times(u1,-1,min_u1); //make here -u1, because you need to flip the u1 vector, to find p1 (and patch angle) correctly.
    p1= find_directing_patchvector(psl, ipart,  min_u1);
    p2= find_directing_patchvector(psl, jpart,  u1);
    
    calc_angles( u1,p1, p2,costheta_i,costheta_j,costheta_ij, calc_ij); // this function gives back theta_i, theta_j and theta_ij (if calc_ij==1)

    return;
}

void calc_angles( vector u1,vector p1,vector p2, double *ptrcostheta_i, double *ptrcostheta_j,double *ptrcostheta_ij, int calc_theta_ij){
    /* calclulates the angles between a interparticle vector u1 and two patch vectors p1 and p2*/
    double costheta_i,costheta_j,costheta_ij,rad_to_degree=180./PI;
    
    
    costheta_i=-vector_inp(u1,p1); 
    costheta_j=vector_inp(u1,p2); 
    // gprint(costheta_i);
    // gprint(costheta_j);
    // vprint(u1);
    // vprint(p1);
    // vprint(p2);


    *ptrcostheta_i=costheta_i;
    *ptrcostheta_j=costheta_j;

    if (costheta_i>=1 || costheta_j>=1 || calc_theta_ij==0){
        // *theta_i=(costheta_i<=1.)? acos(costheta_i)*180./PI:0;
        // *theta_j=(costheta_j<=1.)? acos(costheta_j)*180./PI:0;
        // *theta_ij=0.;
        *ptrcostheta_ij=1;
        
    } 
    else{ 
        vector e1,e2,e3,proj_u1_p1,proj_u1_p2,u2,u3; 

        // projection p1 onto u1, and create (unit vector) e2 the perpendicular projection of p1 
        scalar_times(u1,-costheta_i,proj_u1_p1); 
        vector_minus(p1,proj_u1_p1,u2);
        scalar_divide(u2,sqrt(vector_inp(u2,u2)),e1);
        // projection p2 onto u1, and create (unit vector) e3 the perpendicular projection of p2 
        scalar_times(u1,costheta_j,proj_u1_p2); 
        vector_minus(p2,proj_u1_p2,u3);
        scalar_divide(u3,sqrt(vector_inp(u3,u3)),e2);

        /*calc angle3={0,180}, you need it in angles, because the interpolation needs angles*/
        *ptrcostheta_ij=vector_inp(e1,e2);

        //return the costheta, not the theta
        // *theta_i=acos(costheta_i);
        // *theta_j=acos(costheta_j);
        // *theta_ij=(costheta_ij);
    }

    return;
}

double potential_attractive_bond_energy(Slice *psl, int ipart, int jpart, double s, vector u1){
    /*this fucntion calculates Vc*S and returns this energy
    you already know that s is small enough < rccutoff
    */
    int typei,typej;
    double cositheta,cosjtheta,cosijtheta,S, Eatr=0, Ebond=0;
    double cosdelta_i,cosdelta_j;
    

    typei=psl->pts[ipart].ptype;
    typej=psl->pts[jpart].ptype;

    // PUT THIS IN INIT AT MAXDEFS CHECK
    if ((typei>sys.nparticle_types) || (typej>sys.nparticle_types)){
        dprint(typei);
        dprint(typej);
        dprint(sys.nparticle_types);
        error("potential_attractive_bond_energy: the Ptype is too large");
    }
    
    //the code is ready for isotropic particles: make cosdelta=-1 (180 graden) and s_fixed 1 ; add as many patches/sites you like as a maximum of bonds
    cosdelta_i=site[typei].cosdelta; 
    cosdelta_j=site[typej].cosdelta;
   
    orientation_parameters( psl,  ipart,  jpart, u1, &cositheta, &cosjtheta, &cosijtheta);

    /*consider the treshold of the patches */
    // if ((sys.particletype[typei].nsites==3) && (sys.particletype[typej].nsites==3)){
    //     S=0; // no interaction between two tripatch particles
    //         //WARNING THIS IS SPECIFIC FOR ACTIVE NETWORK
    // }else 
    if ((cositheta<cosdelta_i) || (cosjtheta<cosdelta_j)){
        S=0;
        return 0;
    }
    else{
        S=S_value(cositheta,cosjtheta, cosijtheta, typei, typej);
        // gprint(S);
    }

    

    /*The attractive part of Simon Potential*/
    /*here use an -A/xi exp(-r^2/xi) function where A and xi are dependent on the r_wetting and surface_charge*/
    Eatr = potential_attractive_energy_sdist(s);

    Ebond = S* (Eatr);
    // if (Ebond<0){
    //     dprint(ipart);
    //     dprint(jpart);
    //     gprint(S);
    //     gprint(Eatr);
    //     vprint(psl->pts[ipart].r);
    //     vprint(psl->pts[jpart].r);
    //     gprint(cositheta);
    //     gprint(cosjtheta);
    // }

    if(Ebond!=0.){
        if (Ebond/Ebond!=1 ){
            gprint(S);
            gprint(cositheta);
            gprint(cosjtheta);
            error("Ebond is nan");

        }
    }
    return Ebond;
}

void energycheck_nn(Slice *psl){
    int i, nearest_neighbor_old=sys.nearest_neighbor;

    double en_wo,  en;
    en=total_energy(psl);

    sys.nearest_neighbor=0;
    en_wo=total_energy(psl);
    sys.nearest_neighbor=1;
    // printf("end\n\n");
    if(en-en_wo>=1e-6){
        printf("WARNING energy %lf not equal to energy with neighborlist %lf ", en_wo, en );
        printf("with a difference of     %lf\n", en_wo-en);

        // en=potential_energy(psl);
        if(fabs(en-en_wo)>=0.9  ){
            error("The difference is larger than 1.0 kT; quitting using the nearestneighbor list\n");
        }
    }
    sys.nearest_neighbor=nearest_neighbor_old;
    return;
}

double gravitational_energy_ipart(Slice *psl, int ipart){
    double gravitational_energy=0.0, zinv,zinv2,zinv6,zinv12,psiz;
    double fg,zcut,b_zc,radius;
   
    Pts *psi;
    psi = &psl->pts[ipart];

    
    zcut=sys.particletype[psi->ptype].zcut;
    fg=sys.particletype[psi->ptype].fg;
    b_zc=sys.particletype[psi->ptype].b_zc;
    
    radius=sys.particletype[psi->ptype].radius;
    
    psiz =psi->r.z-radius;


    if(psiz<0){
        psiz+=sys.boxl.z;
    }
    if(psiz>=zcut){ 
        gravitational_energy = fg*psiz - b_zc; /*V_g = Fg*z - b*/
    }
    else{
        zinv = 1./psiz;
        zinv2 = zinv*zinv;
        zinv6 = zinv2*zinv2*zinv2;
        zinv12 = zinv6*zinv6;
        gravitational_energy = 4.*pot.epsilongravLJ*(zinv12-zinv6+1./4.);
    }
    // gprint(zinv12);
    // gprint(zinv6);
        

    return gravitational_energy;


}

void particle_pair_energy(Slice *psl, int ipart, int jpart, double *Erep, double *Ebond, double *s){
    // calculates the terms of the pair potential between ipart and jpart. it returns Erep, Ebond and s. 
    // which may be used to save bond information.

    double r;
    vector rij,u1;

    particles_distance_length_vector(psl,ipart,jpart,s,&r,&rij);

    if((*s<pot.s_cutoff )&& (*s>pot.s_overlap)) { 
        // printf("particles are making are close .. ");

        //calculate repulsive energy(=isotropic,i.e. not dep on patch angles)
        *Erep=potential_repulsive_energy_sdist(*s);
        // gprint(*Erep);
        if(isnan(*Erep)!=0){
            gprint(*Erep);
            error("is nan");
        }
        /*this calculates: Vc(r)*S */
        scalar_divide(rij,r,u1);
        *Ebond=potential_attractive_bond_energy(psl,ipart,jpart,*s,u1);
        // gprint(*Ebond);
        
        return;
    }
    else if(*s<=pot.s_overlap){
        *Erep=1e6;
        *Ebond=0.;
        return ;
    }
    else{
        *Ebond=0.;
        *Erep=0;
    }    
    return;
}

double total_Epair_neighborlist(Slice *psl){
     
    /* the potential energy using the neighbor list
      the neighbor list is a structure that works with 26 images */
 
    Pts *psi,*psj;

    double s;
    double Erep=0,Ebond=0;
    double Erep_tot=0,Ebond_tot=0 ,Epair=0;

    int i,k,a,j,count,*b,tot_bonds=0;

    // printf("*** with negihborlist***\n");

    
    //calculate patchy Simons potential interaction between the neighboring particles 
    /* use neighborlist from here on*/
    count=0;

    for(a = 0;a<list.nimages;a++) {
        i = list.nli[a].ipart;
        k = list.nli[a].image_id;

        if ((k<0)|| (k>=NSHIFTS)) error("force:k out of range\n");   
        // psi = &psl->pts[i];
        // vector_add(psi->r,list.trans[k],ri) ; // i think this one is neccessary when using the cells

        for(b=list.nli[a].first; b<list.nli[a+1].first; b++){
            j = *b;                     /* b is a pointer to the next particle number in the image */
            // psj =&psl->pts[j];  

            particle_pair_energy(psl,i,j,&Erep, &Ebond, &s);

            Erep_tot += Erep; 
            Ebond_tot += Ebond ;

            if((Ebond)<pot.bond_cutoffE){
                // n eighborlist has only unique combinations
                 // printf("there is a bond: with Vpair=%lf\n", Ebond+Erep);
                add_bond_information(psl,i,j,Erep,Ebond,s);                    
                add_bond_information(psl,j,i,Erep,Ebond,s);   
                tot_bonds++;
            }
            
            count++;
        }
    }
    if (count!= list.nneighbors) error("nlist corrupted in energy"); 
    // printf("tot bonds are %d \n\n",tot_bonds);
    psl->nbonds=tot_bonds;
    Epair = Ebond_tot + Erep_tot ;

    return Epair;
}


double total_Epair_woneighborlist(Slice *psl){
    int ipart,tot_bonds=0;
    double Erep_tot=0,Ebond_tot=0 ,Epair=0;
    double s,Erep=0,Ebond=0;

    // printf("*** total energy without  neighborlist***\n");

    //calculate patchy critical Casimir interaction between all particles
    for(  ipart=0; ipart<psl->nparts; ipart++ )  {
        /*kpart is de teller, bereken de energy van ipart, start counting the sum vanaf jpart*/
        for( int kpart=ipart+1; kpart<sys.npart; kpart++ ) {
      
            particle_pair_energy(psl,ipart,kpart,&Erep, &Ebond, &s);
            
            Erep_tot += Erep; 
            Ebond_tot += Ebond ;

            if(Ebond<pot.bond_cutoffE){
                // printf("there is a bond (%d,%d) with Vpair=%lf\n", ipart,kpart,Ebond+Erep);
                add_bond_information(psl,ipart,kpart,Erep,Ebond,s);
                add_bond_information(psl,kpart,ipart,Erep,Ebond,s);   
                tot_bonds++;   
            }
        }

        // tot_bonds+=psl->pts[ipart].nbonds;
    }
    // printf("particle energies done\n");
    psl->nbonds=tot_bonds; // bonds are counted twice, so divide by 2

    Epair = Ebond_tot + Erep_tot ;

    return Epair;
}
