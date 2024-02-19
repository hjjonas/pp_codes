#include "path.h"

/*------------------LOCAL FUNCTIONS------------------------------------------*/
void wall_interaction_force_ipart(Slice *, int );
void activity_force_torque(Slice *, int );
void integrate_time(Slice * , double , double , int *);
void gravitational_force_ipart(Slice *, int );
void single_bond_force(Slice *, int , int);
void calculate_total_forces(Slice *);
void calculate_forces_neighborlist(Slice *);
void gravitational_force_ipart(Slice *, int );
void derivative_check(Slice *);
void propagate_bd(Slice *) ;
/*---------------------------------------------------------------------------*/


void bmdcycle(Slice *psl) {
    /*sys.sim_type==0 Brownina Dynamics*/
    int istep, update=0, ipart;
    double dr2;

    // Perform overdamped langevin dynamics for langevin.step steps
    for(int istep=0;istep<langevin.step;istep++){
        propagate_bd(psl);        
    }

    // If bond tracking is enabled, perform bond tracking for the current slice
    // if(analysis.bond_tracking==1){
    //     bond_tracking(psl);
        // printBondTimeInfoArray(bondtimeinfoarray )
    // }
    if (analysis.print_trajectory) printing_trajectory(&slice[0]); // try to do it every 1 second through langevin.ninter*timestep

}

void propagate_bd(Slice *psl)  {

    int update=0;
    
    for( int inter=0; inter<langevin.ninter; inter++) {
        /* calculate the forces on all particles, and integrate in time*/
        calculate_total_forces(psl); 
        integrate_time(psl,  langevin.timestep, langevin.dtD, &update);
        
        if ((sys.nearest_neighbor==1 )&& (update)) {
            // printf("Updating the neighborlist\n");
            update_nnlist(psl);
            update =0;
        }
    }

    // cumulative time 
    psl->c_time+=(langevin.timestep*langevin.ninter);

    return;
}

void integrate_time(Slice *psl, double timestep, double timestep_randomD, int *nn_update){
    // timestep is input parameter, and timestep_randomD is sqrt(2*timestep/beta) (see init.c)
    Pts *psi; 
    Particletype *ptypei;

    vector xi,f,t,u,fa;
    quaternion qu1,qu2,qu3,qprime,qua;
    tensor rotmat;
    quattensor Bmat;

    double dt,lambdaq,fnorm, dr2;


    for( int ipart=0; ipart<psl->nparts; ipart++) {
        psi = &psl->pts[ipart];
        ptypei=&sys.particletype[psi->ptype];

        //TRANSLATION part due to force
        //DT*F*dt*Beta = muT*F*dT see eqn 9 first term of Ilie2015
        scalar_times(psi->f,timestep,f);
        scalar_times(f,langevin.mobilityT,f);
        vector_add(psi->r,f,psi->r);
        if(sys.nearest_neighbor){ vector_add(psi->dr,f,psi->dr); }

        //random translation
        xi=RandomBrownianVector(timestep_randomD); // random fluctuation uses this langevin.dtD=sqrt(2*timestep/beta)
        scalar_times(xi,langevin.sqrtmobilityT,xi); 
        vector_add(psi->r,xi,psi->r);
        if(sys.nearest_neighbor){ vector_add(psi->dr,xi,psi->dr); }
        


        if(sys.nearest_neighbor==1){ 
            dr2 = vector_inp(psi->dr , psi->dr);
            if ((dr2 > list.cutoff_buffer2)) {
                *nn_update=1;
            }
        }
        

        //ROTATION
        //get matrices for q(t)
        //rotmat=A (eqn 12 Ilie2015), Bmat= matrix for quaternion rotation: Baalpha
        rotmat = getrotmatrix(psi->q); 
        Bmat = getquatmatrix(psi->q);

        //rotational part due to force. eqn 13 first term is qu1
        //Baalpha * (muR) * rotmatA * torque * delt = qu1
        // vprint(psi->t);
        scalar_times(psi->t,timestep,t);
        //do not forget, multiply with the inverse rotation matrix to convert to body-fixed torque
        matrixT_x_vector(rotmat,t,u);
        /*unit of mobility_R, rad/s? anwser: I think it is in units of sin(a/2) where a is the rotation angle.

        first, vector(torque) x scalar -> vector(3x1) , matrix(3x3) x vector(3x1) -> vector(u), 
        vector u (3x1)x scalar mobilityR (1x1)-> vector(3x1)  , matrix Bmat (4x3) x vector(3x1) -> quaterion (4x1)
        it leads to a quaternion q.0 = cos(a/2) q.1 =sin(a/2)i, q.2 = sin(a/2)j, q.3=sin(a/2)k where (i,j,k) is the rotating axis.
        as you start with a vector --> (i,j,k) and you multiply it with
            */
        scalar_times(u,langevin.mobilityR,u);
        /* ok here you multiply u (torque x timestep x drag) */
        quatmatrix_x_vec(Bmat,u,qu1);
        
        //Baalpha * (muR) * xi = qu2
        xi=RandomBrownianVector(timestep_randomD);
        scalar_times(xi,langevin.sqrtmobilityR,xi);
        //Bmat is the matrix representation of the quaternion of the particle at time t
        quatmatrix_x_vec(Bmat,xi,qu2);

        //qprime = qu1+qu2+q(t) for the first time
        quat_add(qu1,qu2,qprime);
        quat_add(qprime,psi->q,qprime);

        //find lambdaq, I guess it is th smallest...check which one keeps q(t+delt)^2=1
        //woohoo it works for the smaller lambdaq, maybe also simply for the bigger...
        lambdaq=langrange_multiplier_quat(qprime, psi->q);
        if(lambdaq>1e20) {
            //langrange multiplier did not work, simply renormalize quaternion
            scdivide_quat(qprime,sqrt(quat_inp(qprime,qprime)),psi->q);
        }
        else {
            //lambdaq*q(t) = qu3
            sctimes_quat(psi->q,lambdaq,qu3);
            //q(t+delt) = qprime+qu3
            quat_add(qprime,qu3,psi->q);
        }
    }
    // update the "current" patch vectors after you updated the quaternions. 
    //This update is specifically necessary when you want to calculate the new forces and torques 
    update_patch_vectors(psl);

    return;
}

/*__________________FORCE EXPRESSIONS_________________*/
void calculate_total_forces(Slice *psl) {
    
    Pts *psi;
    Particletype *ptypei;
    int ipart,jpart;

    //initialize all forces to zero and update patch vectors
    for( ipart=0; ipart<psl->nparts; ipart++ )  {
        psi = &psl->pts[ipart];
        psi->f = nulvec;
        psi->t = nulvec;
         // put all nbonds to zero to restart the tracking of the bonds
        psi->nbonds =0;
    }

    //calculate forces based on force definition is simonsparameters.c
    if (sys.nearest_neighbor==1){
        calculate_forces_neighborlist(psl);
    }
    else{
        for( ipart=0; ipart<psl->nparts; ipart++ ) {
            for( jpart=ipart+1; jpart<psl->nparts; jpart++ ) {
                single_bond_force(psl,ipart,jpart);
            }
        }
    }

        // activity and gravity
    for( ipart=0; ipart<psl->nparts; ipart++ )  {
        psi = &psl->pts[ipart];
        ptypei=&sys.particletype[psi->ptype];

        if(ptypei->activity==1){
            // deterministic self-propulsion for translation
            // this part of the force is 
            activity_force_torque(psl,ipart);
        }
        if(sys.gravity>0){
            gravitational_force_ipart(psl,ipart);

            if (pot.wall_int!=0){
                wall_interaction_force_ipart(psl,  ipart);
            }
        }     
    }

    return ;
}

void activity_force_torque(Slice *psl, int ipart){
    /* add here the activity force on the particle ipart*/
    
    Pts *psi;
    psi = &psl->pts[ipart];
    vector fa,Fa;
    double rad_to_degree=180./PI;
    int ptype=psi->ptype;
    
    tensor rotmat = getrotmatrix(psi->q); 
    matrix_x_vector(rotmat,sys.particletype[ptype].e_A,Fa); // e_A unit vector; Fa the direction of the active force
    
    // I added this if statement, to enforce the alignment with the wall but without activity.
    // this is a small hack to simulate the passive system with wall alignment, but without activity in BMD
    // its not very beautifully written, but just a quick hack. 
    // Else I have to add more variables that checks if particles need alignment that are not active. 
    if (sys.particletype[ptype].F_A>0.0001){
        scalar_times(Fa,sys.particletype[ptype].F_A,fa); // multiply with size of F0
        vector_add(psi->f,fa,psi->f); // add to force on psi
        // printf("adding active translation");
    }
    // even if sys.particletype[ptype].F_A<0.0001 then still you execute the alignment below.


    // I want to add alignment of the active force wrt to the wall to prevent flying particles against gravity
    // preform this via the torque
    // its  parabola with U=anlge^2 wrt the xy-plane (or its normal: the z-axis); this torque so not dependent on the magnitude of Fa, anlge in rad
    if (sys.gravity>0){
        double cosangle,angle,magnitude;
        vector facrossz=nulvec,torque=nulvec,proj_facrossz,Fa_proj;

        // the angle via asin
        angle=fabs(asin(Fa.z));

        magnitude=500.*angle;
    
        if( magnitude>0){
            // printf("adding active alignment");
            Fa_proj=Fa;
            Fa_proj.z=0;
       
            normvec(Fa_proj,Fa_proj); //Fa_proj is the projection of Fa onto the xy plane
            vector_cross(Fa,Fa_proj,facrossz); // note that its Fa x Fa_proj (right hand rule)
            normvec(facrossz,facrossz); // facrossz is now the direction of the toruqe vector, 
            
            scalar_plustimes(facrossz,magnitude,torque); // multiply the direction of the torque with the magnitude
            vector_add(torque,psi->t,psi->t);
            
            // printf(" >>>   \n" );
            // printf(" arcsin(e,z) = %lf \n",asin(Fa.z) );
            // printf(" length Fa  = %lf  (  %lf  %lf  %lf )\n",sqrt(vector_inp(Fa,Fa)), Fa.x, Fa.y, Fa.z);

            // printf(" angle = %lf \n",angle );
            // printf(" angle - arcsin = %lf \n",angle - fabs(asin(Fa.z)) );
            if (isnan(torque.x)){
                vprint(Fa);
                vprint(Fa_proj);
                vprint(facrossz);
                vprint(torque);
                vprint(psi->t);

                gprint(angle);
                gprint(magnitude);
                error("stop; torque nan?");
            }
        }
    }
    
    return;
}


void gravitational_force_ipart(Slice *psl, int ipart){
    double zinv,zinv2,zinv6,zinv7,zinv12,zinv13,psiz,fgrav;
    double fg,zcut,b_zc,radius;
    vector z=nulvec;

    Pts *psi;
    psi = &psl->pts[ipart];

    zcut=sys.particletype[psi->ptype].zcut;
    fg=sys.particletype[psi->ptype].fg;
    b_zc=sys.particletype[psi->ptype].b_zc;
    radius=sys.particletype[psi->ptype].radius;
    psiz =psi->r.z-radius;
   
    /*the z unitvector; used for the giving gravity a direction*/
    z.z=1.; 
    if(sys.gravity>0){

        if(psiz<0){
            psiz+=sys.boxl.z;
        }
        if(psiz>=zcut){ 
            fgrav = fg; /*V'_g(z) = -Fg  minus goes later*/
        }
        else{
            zinv = 1./psiz;
            zinv2 = zinv*zinv;
            zinv6 = zinv2*zinv2*zinv2;
            zinv7 = zinv6*zinv;
            zinv12 = zinv6*zinv6;
            zinv13 = zinv12*zinv;
            fgrav = 4.*pot.epsilongravLJ*(-12.*zinv13+6.*zinv7); /* LJ12-6' = 4*epsilonLJ*(-12/r13+6/r7)*/
        }
        
    }

    scalar_mintimes(z,fgrav,psl->pts[ipart].f);
    return ;
}

   
void single_bond_force(Slice *psl, int ipart, int jpart){
    /*the forces acting on the particle: brownian motion and bond forces (no gravitational force)
    Allen paper: 10.1080/00268970601075238 Expressions for forces and torques in molecular simulations using rigid bodies, allen & germano */

    Pts *psi, *psj;   
    vector rij;
    double r,s, rinv;
    int typei, typej;
    
    psi = &psl->pts[ipart];
    psj = &psl->pts[jpart];
    typei=psi->ptype;
    typej=psj->ptype;


    particles_distance_length_vector(psl,ipart,jpart,&s,&r,&rij); 
    // printf(" via bond distance fuctinn s=%lf , r=%lf ,rij  = (%lf  %lf  %lf)\n",s,r,rij.z,rij.y,rij.z);
    // vprint(psl->pts[ipart].r);
    // vprint(psl->pts[jpart].r);



    if(s<=pot.s_overlap ){
        printf("s = %lf <= %lf \n", s,pot.s_overlap);

        printf("r = %lf \n", r);
        vprint(rij);
        dprint(ipart);
        vprint(psi->r);
        dprint(jpart);
        vprint(psj->r);
        error("Overlapping particles!!!");
    } 


    // printf("the calculation the force between particles %d and %d starts \n",ipart,jpart );
    if(s>pot.s_overlap && s<pot.s_cutoff ){ 
        double fmag,S;
        double cositheta,cosjtheta,cosdelta_i,cosdelta_j,cosijtheta;
        vector rnorm,p1,p2,min_rnorm;
        
        // repulsive radial force
        fmag=fabs(bond_repulsive_force(s));/*F=-dVrep/ds >= 0*/

        scalar_divide(rij,r,rnorm); 
        // vprint(rnorm);

        /*add the repulsive force (=a positive number) positively to i and negatively to j (as rij points toward from j to i)*/
        scalar_plustimes(rnorm,fmag,psi->f);
        scalar_mintimes(rnorm,fmag,psj->f);
        /*repulsive radial force done*/
    
           

        if ((sys.particletype[typei].nsites==3) && (sys.particletype[typej].nsites==3)){
            S=0; // no interaction between two tripatch particles
            //WARNING THIS IS SPECIFIC FOR ACTIVE NETWORK
        }
        else{
            //the part written below, upto calc_angles, is a copy of orientation_parameters() in energy.c
            //find smallest combination of the patch vector angles
            orientation_parameters( psl,  ipart,  jpart, rnorm, &cositheta, &cosjtheta, &cosijtheta);
            S = S_value( cositheta,  cosjtheta,cosijtheta  , typei,  typej);

        }
        
        
        if (S>0){
            double Umag,fmagP,fmagrinv,dSdcostheta;
            vector rcrosspi,rcrosspj,piperpr,pjperpr;

            /*if you are here. you found the patches that are connected*/
            // add_bond_information(psl,ipart,jpart,0,0,0);
            // add_bond_information(psl,jpart,ipart,0,0,0);

            // /* attractive energy and force of the patches wihtout angular correction */
            Umag = fabs(potential_attractive_energy_sdist(s)) ;  // Umag is a negative number; do absolute value
            fmag = fabs(bond_attractive_force(s)); // Fc <=0 , thus fabs() 

            // the magnitude of the force on the patches. take absolute to get magnitude
            fmagP = S*fmag;

            /*the effective attractive force(fmagP). The particle attract, so mintimes on i and plustimes on j */

            //the radial attractive force  along rnorm
            scalar_mintimes(rnorm,fmagP,psi->f);
            scalar_plustimes(rnorm,fmagP,psj->f);

            // printf("patch radial interaction \n");
            // vector cross is right hand rule, wijsv. = a , middel v. = b, duim: a x b, 
            /*the torque = p x rnorm , so here rcrosspi = -torque (due to direction of rnorm)*/
            scalar_times(rnorm,-1.,min_rnorm); //make here -rnorm, because you need to flip the rnorm vector, to find p1 (and patch angle) correctly.
            p1= find_directing_patchvector(psl, ipart,  min_rnorm);
            p2= find_directing_patchvector(psl, jpart,  rnorm);

            vector_cross(rnorm,p1,rcrosspi);
            vector_cross(rnorm,p2,rcrosspj);

            /*the direction of the force that bring patch back to center,ie. aligning with rnorm*/
            /*rnorm points form j to i*/
            /*pi is projected onto pi*/
            vector_cross(rnorm,rcrosspi,piperpr);
            vector_cross(rnorm,rcrosspj,pjperpr)

            //"sliding" radial force, the one perperdicular to rnorm see eqn 4 of the Allen paper
            //CALCULATE DEL U / DEL COSTHETA_i derivative to i
            dSdcostheta = dS_dcostheta(psl,cositheta, ipart,cosjtheta,jpart);
            fmag = fabs((Umag)*dSdcostheta);
            rinv=1./r; 
            fmagrinv=fmag*rinv;

            scalar_plustimes(piperpr,fmagrinv,psi->f);
            scalar_mintimes(piperpr,fmagrinv,psj->f);
            // vprint(psi->f);
            // vprint(psj->f);


            /*rcrosspi is -torque and fmag is -F, therefore plustimes */
            // TORQUE on particle i
            // see after eq. 4 of Allen Paper for torque expression
            scalar_plustimes(rcrosspi,fmag,psi->t);
            
            // DEL U / DEL COSTHETA_j derivative to j
            dSdcostheta  = dS_dcostheta(psl,cosjtheta,jpart,cositheta, ipart);
            fmag = fabs((Umag)*dSdcostheta);
            fmagrinv=fmag*rinv;
            scalar_mintimes(pjperpr,fmagrinv,psi->f);
            scalar_plustimes(pjperpr,fmagrinv,psj->f);

            // TORQUE on particle j
            scalar_mintimes(rcrosspj,fmag,psj->t);


        }

    }


    return;
}
void calculate_forces_neighborlist(Slice *psl){
    /*uses the neighbor list to loop over the particles and calculates the single bond force
    the gravitational force and patch vectors are already calculated at this point*/

    // Pts *psi,*psj;
    // vector rij,ri;  
    int i,k,a,j,*b;

    list.count=0;
    for(a = 0;a<list.nimages;a++) {
        i = list.nli[a].ipart;
        k = list.nli[a].image_id;
        if ((k<0)|| (k>=NSHIFTS)) error("force:k out of range\n");   
        // psi = &psl->pts[i];
        /*vector_add(psi->r,list.trans[k],ri) ; not neccessary , I use option (1) in neighborlist not using the cells*/

        for(b=list.nli[a].first; b<list.nli[a+1].first; b++){
            j = *b;
            // psj =&psl->pts[j];  not neccessary for the current funtion of single_bond_force


            // give the rij to the single_bond_force ; 
            // single_bond_force( psl,  i,  j,rij); //  with rij vector
            single_bond_force( psl,  i,  j);

            list.count++;
        }
    }
    if (list.count!= list.nneighbors) error("nlist corrupted in force"); 

    return;
}
double F_harmonic_oscillator(double k, double x, double x_shift){
    return -1*k*fabs(x-x_shift);
}
void wall_interaction_force_ipart(Slice *psl, int ipart){
    /* if the wall has an interaction with the patches
    there is a harmonic oscillator that governs the interaction strength with its minimum below 0*/
    Pts *psi;
    psi=&psl->pts[ipart];
    int isite;
    vector z_axis=nulvec,x_axis=nulvec,y_axis=nulvec,p1;
    double Fwall_tmp=0,Fwall,z,angle_max=site[psi->ptype].cosdelta;
    double zcut=sys.particletype[psi->ptype].zcut+sys.particletype[psi->ptype].radius*1.1,angle,cosangle,S,F_height;
    z_axis.z=1;
    x_axis.x=1;
    y_axis.y=1;

    vector rcrosspi,piperpr;
    double F_harm, E_harm, fmag, Umag,fmagP,dSdcostheta,fmagrinv;

   
    // only if the particle is in the vicinity of the wall: 
    if (psi->r.z<=zcut & psi->r.z>0){
        // find thenpatch that is pointing downwards. take the innerproduct with patch vector if 
        for(isite=0;isite<sys.particletype[psi->ptype].nsites;isite++){
            p1=psi->patchvector[isite];
            if (p1.z>0){
                continue;
            }
            cosangle=-vector_inp(p1,z_axis); 
            // angle=acos(cosangle);
            if (cosangle<angle_max){
                continue;
            }
            
            // F_height=gravitational_force_ipart(psl, ipart); // gravitational energy is zero at the wall
            // z=fabs(p1.z-zcut);

            Umag=fabs(harmonic_oscillator(pot.wall_int,p1.z,zcut));
            fmag=fabs(F_harmonic_oscillator(pot.wall_int,p1.z,zcut));
            
            S= Saccent( cosangle, psi->ptype);

            fmagP = S*fmag;
            
            scalar_mintimes(z_axis,fmagP,psi->f);

            /*the torque = p x rnorm , so here rcrosspi = -torque (due to direction of rnorm)*/
            vector_cross(z_axis,p1,rcrosspi);

            /*the direction of the force that bring patch back to center*/
            /*rnorm points form j to i*/
            /*pi is projected onto pi*/
            vector_cross(z_axis,rcrosspi,piperpr);
            
            //"sliding" force, the one perperdicular to rnorm
            //CALCULATE DEL U / DEL COSTHETA_i derivative to i
            
            dSdcostheta = dSaccent_dcosangle(cosangle,psi->ptype);
            fmag = fabs((Umag)*dSdcostheta);
            fmagrinv=fmag/z;

            scalar_plustimes(piperpr,fmagrinv,psi->f);

            /*rcrosspi is -torque and fmag is -F, therefore plustimes */
            scalar_plustimes(rcrosspi,fmag,psi->t);
        }
    }
    return  ;
}
void forcecheck_nn(Slice *psl){
    /* testing the force difference with and without using neighborlist*/
    vector forcedif=nulvec, tracker=nulvec;
    int i, ipart, jpart;
    vector af[psl->nparts], bf[psl->nparts];
    vector ar[psl->nparts], br[psl->nparts];
    double df2;

    // printf("checking the force of nn is good \n");
    for(i=0;i<psl->nparts;i++){
        bf[i]=nulvec;
        af[i]=nulvec;
        br[i]=nulvec;
        ar[i]=nulvec;

        psl->pts[i].f=nulvec;
    }

    //calculate forces based on force definition is simonsparameters.c
    for( ipart=0; ipart<psl->nparts; ipart++ ) {
        for( jpart=ipart+1; jpart<psl->nparts; jpart++ ) {
            //just give here nulvec as vector, it rij will be calculated if sys.nearest_neighbor!=1
            // single_bond_force(psl,ipart,jpart, nulvec);
            single_bond_force(psl,ipart,jpart);
        }
    }
    
    for(i=0;i<psl->nparts;i++){
        bf[i]=psl->pts[i].f;
        br[i]=psl->pts[i].r;
        psl->pts[i].f=nulvec;
    }

    calculate_forces_neighborlist(psl);
    for(i=0;i<psl->nparts;i++){
        af[i]=psl->pts[i].f;
        ar[i]=psl->pts[i].r;
    }

    for(i=0;i<psl->nparts;i++){
        vector_minus(af[i], bf[i], forcedif );
        vector_add(forcedif, tracker, tracker)
    }
    // vprint(tracker);
    df2 = vector_inp(tracker,tracker);
    if(fabs(df2)>0.9){
        printf("The df2 is %.12lf \n", vector_inp(tracker,tracker));
        error("force of neighborlist not good");
    }
    // else{
    //     // printf("Force is good, difference df2 is %lf\n", df2);
    // }
    return;
}

/*__________________INITIALIZA BMD_________________*/
void setup_BMD(Slice *psl){
    //do here the force check, you need s_min for that
    printf("\nChecking the derivatives/forces...");
    derivative_check(psl);
    printf("done\n");
   
   /*langevin/brownian md parameters*/
    printf("\nDefining the langevin timestep and mobility parameters... ");

    // so the input used to read mobility and beta, but its odd. 
    // I changed it to reading in diffusion 13-okt, so mobility is calculated inside the code via beta
    // https://en.wikipedia.org/wiki/Einstein_relation_(kinetic_theory)
    // D = mu*kT = mu/Beta ; thus mu = D*beta
    if (langevin.mobilityT <1e-5 || langevin.mobilityR <1e-5){
        gprint(langevin.mobilityT);
        gprint(langevin.mobilityR);
        error(" the rotational diffusion and/or translational diffusion/mobiilty is close to zero? ");
    }
    langevin.diffusionT= langevin.mobilityT/sys.beta;
    langevin.diffusionR= langevin.mobilityR/sys.beta;
    
    langevin.sqrtmobilityT = sqrt(langevin.mobilityT); //eqn 5,  
    langevin.sqrtmobilityR = sqrt(langevin.mobilityR); // eqn 4 Ilie2015

    langevin.dtD=sqrt(2.0*langevin.timestep*psl->temp);
    //langevin.dtBeta=langevin.timestep*psl->beta; // oeeehhh I see that beta is put in here... 

    
    // depending on the input you give, you determine the ncycle1, langevin.step
    if (sys.ncycle1==0){
        langevin.step=(int)ceil(langevin.print_time/(langevin.timestep*langevin.ninter));
        sys.ncycle1=langevin.total_time/(langevin.step*langevin.timestep*langevin.ninter*sys.ncycle2);
    }
    else{
        langevin.step=100;
    }

    if(init_conf.restart){
        restart_read_time(psl);
    }

    if (sys.ncycle1==0 || sys.ncycle2==0 || langevin.step==0 || langevin.ninter==0){
        dprint(sys.ncycle1);
        dprint(sys.ncycle2);
        dprint(langevin.step);
        dprint(langevin.ninter);
        error(" one of these is zero. stop");
    }
    
    if (init_conf.mc_warmup>0){
        sprintf(mc_single.name,"single particle moves");
        setup_mc_move(&mc_single);

    }

    printf("done\n");

        /* create the neighbor list here;*/
    if(sys.nearest_neighbor==1){
        printf("   step 1) Setting up the neighbor list\n");
        setup_nnlist();
        printf("   step2) filling up the nnlist  \n");
        update_nnlist(psl);
    }

    return;
}

void derivative_check(Slice *psl){
    /* this is to compare (and check) the angular part and (radial) of derivative the potential
    The derivative is numerically approximated at various delta_cos and the difference with the exact solution.

    I check the derivative at max_iter_q=100 points (0.25 degrees from eachother)
    so at theta=0,1,2,3,4,..., 25 degrees
    I calculate the energy at theta, theta+delta_cos etc where delta_theta= E-6, E-5, E-4, E-3, E-2, E-1
    Calculate the force numerically by evaluating the energy difference F = -dE/dr 
    note thet dr = 0.5PI * delta_theta / 180  

    There are 2 particles,
    rotate one particle, keep the other fixed*/


    double delta_theta=0.1,costheta;
    double dV_NUM, dV_exact,max_angle=25.0, dangle=0.25;
    int max_iter_p = 10, max_iter_q=(int)max_angle/dangle;
    double cos_ijtheta=1;
    
    int p,p0=0,q, old_nsites,gravity_psl,npart_psl, ipart, pt;
    vector old_site, rnorm,p1,p2,rij;
    quaternion dq,dqrot;
    double theta_check,y_prev,dcostheta, Fnum, Fexact,E0,E1,theta=0,r,cos1theta,cos2theta;
    double Eatr, S1,S0;

    Pts *psi, *psj;
    Particletype old_pts_types[PTYPES];

    
    //create a temporary slice for the 2 particles to perform the force check
    Slice *copyslice = malloc(sizeof(Slice));
    memcpy(copyslice, psl, sizeof(Slice));

    // save old settings
    gravity_psl=sys.gravity;
    npart_psl=psl->nparts;

    //turn gravity off and set npart to 2
    sys.gravity=0;
    psl->nparts=2;
    printf("\n");
    //loop over the particle types
    for (pt=0; pt<sys.nparticle_types;pt++){
        // if (site[pt].S_fixed>0){
        //     printf("***constant sites *** \n");
        // }
        // printf("*** site type %d *** \n",pt);
        // printf("   Theta_0  | delta Theta | check Theta |           F num     (S1,S0)           |  F exac  |  Fexact-Fnum   \n");

    
        //save nsites, and first site of particle type
        old_nsites=sys.particletype[pt].nsites;
        old_site=sys.particletype[pt].site[0];
        

        //put the nsites to 1 and site as (0 1 0)
        sys.particletype[pt].nsites=2;
        sys.particletype[pt].site[0]=nulvec;
        sys.particletype[pt].site[0].y=1;

        y_prev=-0.48*sys.boxl.y;
       // setup the 2 particles, which both have 1 patch
        for( ipart=0; ipart<psl->nparts; ipart++) {
            psi=&copyslice->pts[ipart]; 
            psi->ptype=pt;

            psi->r=nulvec;
            psi->r.y = y_prev + (sys.particletype[psi->ptype].radius+pot.s_min) ;                    
            y_prev = psi->r.y+sys.particletype[psi->ptype].radius;
            
            //this makes the patches in opposite directions
            if(ipart==0) {
                psi->q.q0 = 1.0;
                psi->q.q1 = 0.0;
                psi->q.q2 = 0.0;
                psi->q.q3 = 0.0;
            }
            if(ipart==1) {
                psi->q.q0 = 0.0;
                psi->q.q1 = 1.0;
                psi->q.q2 = 0.0;
                psi->q.q3 = 0.0;
            }
        }


        vector_minus(copyslice->pts[0].r,copyslice->pts[1].r,rij);
        r = sqrt(vector_inp(rij,rij));
        scalar_divide(rij,r,rnorm);

        // first loop changes the cos / composition
        theta=0;
        Eatr = potential_attractive_energy_sdist(pot.s_min);
        for(q=0; q<max_iter_q;q++){
            
            // Calculate the energy at delta_theta=0
            update_patch_vectors(copyslice);

            p1=copyslice->pts[0].patchvector[0];
            p2=copyslice->pts[1].patchvector[0];

            cos1theta=-vector_inp(rnorm,p1);
            cos2theta=vector_inp(rnorm,p2);

            costheta=cos(theta/180.*PI);

            
            S0=S_value(cos1theta,cos2theta, cos_ijtheta, pt,pt);
            E0=S0*Eatr;

            // Fexact=dSaccent_dcosangle(costheta, pt);
            Fexact=dS_dcostheta(copyslice,cos1theta, 0,cos2theta,1)*dcostheta;
            // Fexact = dSaccent_dcosangle(cos1theta,pt)*Saccent(cos2theta,pt);

           
            delta_theta=0.1;
            //second loop changes delta_theta= 0.001, 0.01, 0.1 etc
            //I put p0 and max_iter_p to 6 and 7, because I only want to check at delta theta=1E-6
            p0=5;
            max_iter_p=p0+1;
            for(p=0;p<max_iter_p; p++){
                
                //so you start with a conformation at a specific angle 0, 0.25, etc
                //change the angle according to delta_theta , calculate the energy and rotate it back for the next delta_theta
                if (p<p0){ 
                    delta_theta/=10.;
                    continue;
                }
                //rotate particle 0 along x-axis 
                dq=QuaternionXaxis(delta_theta);
                quat_times(dq,copyslice->pts[0].q,dqrot);
                copyslice->pts[0].q = dqrot; 
                update_patch_vectors(copyslice);

                p1=copyslice->pts[0].patchvector[0];
                cos1theta=-vector_inp(rnorm,p1);
                theta_check = cosangle_to_angle(cos1theta);

                // S1=S_value(cos1theta,cos2theta, cos_ijtheta, 0, 0); WRONG the two zero's indicate the particle type
                S1=S_value(cos1theta,cos2theta, cos_ijtheta, pt, pt);

                E1=S1*Eatr;
                
                dcostheta= cos((theta)/180.*PI) - cos((theta+delta_theta)/180.*PI) ;
                Fnum = - (S1-S0);
                // printf("   Theta_0  |   cos1theta  | cos2theta    |   F num     (S1,S0)           |  F exac  |  Fexact-Fnum   \n");

                // printf("    %4.2f    %10.5f  %10.5f  %10.5f (%.5f ,%.5f) %10.5f %10.5f \n", theta,cos1theta,cos2theta, Fnum, S1,S0 ,Fexact,Fexact- Fnum ); //, E0-E1, E0, S0, E1, S1);
                // gprint(S1-S0);
                
                if ((fabs(Fexact- Fnum)>1e-6)){
                    printf("   Theta_0  |   cos1theta  | cos2theta    |   F num     (S1,S0)           |  F exac  |  Fexact-Fnum   \n");
                    printf("    %4.2f    %10.5f  %10.5f  %10.5f (%.5f ,%.5f) %10.5f %10.5f \n", theta,cos1theta,cos2theta, Fnum, S1,S0 ,Fexact,Fexact- Fnum ); //, E0-E1, E0, S0, E1, S1);
                    gprint(cos1theta);
                    gprint(cos2theta);
                    gprint(S1);
                    gprint(S0);
                    gprint(S1-S0);
                    gprint(dcostheta);
                    gprint(dS_dcostheta(copyslice,cos1theta, 0,cos2theta,1));
                    gprint(theta);
                    gprint(delta_theta);

                    gprint(Fnum);
                    error("Fexact- Fnum>1e-3 in the angular force");
                }
                // printf("%20.15f %20.15f %20.15f %20.15f %16.15f %16.15f %20.15f %16.15f %20.15f %16.15f\n", theta, dcostheta , Fnum, Fexact,Fexact- Fnum , E0-E1, E0, S0, E1, S1);
                // printf("    %4.2f       1E-%2d       %20.15f %20.15f %20.15f %20.15f %20.15f \n", theta, p+1 , theta_check, Fnum, E1 ,Fexact,Fexact- Fnum ); //, E0-E1, E0, S0, E1, S1);

                //rotate back particle 0 along x-axis 
                dq=QuaternionXaxis(-1.*delta_theta);
                quat_times(dq,copyslice->pts[0].q,dqrot);
                copyslice->pts[0].q = dqrot; 
                update_patch_vectors(copyslice);

                delta_theta/=10.;
            }
            // printf("________________________________________________________________________________________________________\n");
            dq=QuaternionXaxis(dangle);
            quat_times(dq,copyslice->pts[0].q,dqrot);
            copyslice->pts[0].q = dqrot; 
            update_patch_vectors(copyslice);
            theta=theta+dangle;
            // go from old to new position
        }
        sys.particletype[pt].site[0]=old_site;
        sys.particletype[pt].nsites=old_nsites;
    }

    sys.gravity=gravity_psl;
    psl->nparts=npart_psl;
    /* I checked if original particles are still there, (they are in psl, not in copyslice)*/
    //free the memory of the copyslice
    free(copyslice);

    return;
}

