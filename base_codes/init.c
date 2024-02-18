#include "path.h"



void setup_simulation() {

 
    printf("Setting up the system\n\n");

    printf("initializing the random numbers\n");
    InitializeRandomNumberGenerator(Random_urandom());
    printf("Random_urandom() %d \n", Random_urandom());
    printf("the time(0l) is %ld\n\n", time(0l));

    printf("Allocating memory for the slice\n");
    slice = (Slice *)calloc(1,sizeof(Slice));
    start_slice = (Slice *)calloc(1,sizeof(Slice));
    // // memory_allocation();

    printf("\nReading the input\n");
    read_input(&slice[0]);

    printf("\nInit the model\n");
    init_model(&slice[0]);
    // printf("\nInput read from path.inp:\n");
    // print_input(&slice[0]);


    // if(sys.empty_files==1){
    //     emptyfiles();
    // }

    // //calculates and prints energy, bonding etc
    // printf("\n   calculating energy and bonding of initial configuration:\n");
    // terminate_block(&slice[0]);



    // printf("Setting up the system is done\n");



    return;
}

void read_input(Slice *psl) {

    FILE *fp;
    char *pt,line[MAXLINE];


    sprintf( init_conf.directorypath,"."); // this is used to give the path to e.g. trajectory.xyz 
    if((fp = fopen("path.inp","r"))==NULL) {
        printf("ERROR: could not read input file\n");
        exit(1);
    }

    printf("\nReading input from path.inp\n");
    // Define configuration entries
    // 'd' is an integer, 's' is a string, 'f' is a float, 'p' is are strings, they are just titles to give in input-file structure
    Entry entries[] = {
        //      sys variables
        {"sim_type", &sys.sim_type, 'd'},
        {"ncycle1", &sys.ncycle1, 'd'},
        {"ncycle2", &sys.ncycle2, 'd'},
        {"cluster_MC", &sys.cluster_MC, 'd'},
        {"nearest_neighbor", &sys.nearest_neighbor, 'd'},
        {"switch_method", &sys.switch_method, 'd'},
        {"beta", &psl->beta, 'f'},
        {"boxl.x", &sys.boxl.x, 'f'},
        {"boxl.y", &sys.boxl.y, 'f'},
        {"gravity", &sys.gravity, 'f'},
        // init_conf variables
        {"start_type", &init_conf.start_type, 'd'},
        {"read_path", &init_conf.read_path, 'd'},
        {"directorypath", init_conf.directorypath, 's'},
        {"mc_warmup", &init_conf.mc_warmup, 'd'},
        {"particle_setup", &init_conf.particle_setup, 'd'},
        {"restart", &init_conf.restart, 'd'},
        {"empty_files", &init_conf.empty_files, 'd'},
        {"nchains", &init_conf.nchains, 'd'},
        {"chaingap", &init_conf.chaingap, 'f'},
        //      analysis variables 
        // {"bond_op", &analysis.bond_op, 'd'},
        // {"print_trajectory", &analysis.print_trajectory, 'd'},
        // {"bond_tracking", &analysis.bond_tracking, 'd'},
        // {"rdfanalysis", &analysis.rdfanalysis, 'd'},
        // {"cluster_analysis", &cluster.analysis, 'd'},
        // {"adjacency", &analysis.adjacency, 'd'},
        // {"s_histogram", &analysis.s_distribution, 'd'},
        //  {"xy_print", &analysis.xy_print, 'd'},
        // langevin (BMD) variables
        // {"mobilityT", &langevin.mobilityT, 'f'},
        // {"mobilityR", &langevin.mobilityR, 'f'},
        // {"ninter", &langevin.ninter, 'd'},
        // {"timestep", &langevin.timestep, 'f'},
        // {"print_time", &langevin.print_time, 'f'},
        // {"total_time", &langevin.total_time, 'f'},
        //      potential variables
        {"epsilongravLJ", &pot.epsilongravLJ, 'f'},
        {"dT", &pot.dT, 'f'},
        {"r_wetting", &pot.r_wetting, 'f'},
        {"surface_charge", &pot.surface_charge, 'f'},
        {"wall_int", &pot.wall_int, 'f'},
        {"npart", &sys.npart, 'd'},
        {"s_cutoff", &pot.s_cutoff, 'f'},
        //      'n' denotes a newline
        {"\n", NULL, 'n'}, 
        //      titles in input file, just printed to the screen
        {"SIMULATION\n", "simulation", 'p'},
        {"MC\n", "MC", 'p'},
        {"BMD\n", "BMD", 'p'},
        {"HARMONIC_OSCILLATOR\n", "HARMONIC", 'p'},
        {"FFS\n", "FFS", 'p'},
        {"ANALYSIS\n", "Analysis", 'p'},
        {"GRAVITY\n", "gravity", 'p'},
        {"SWITCH_FUNCTION\n", "switch function", 'p'},
        {"CRITICAL_CASIMIR\n", "critical Casimir potential", 'p'},
        {"SYSTEM\n", "system", 'p'},
        {"K\n", "K", 'p'}
    };

    // read the file and load in the input
    while (fgets(line, MAXLINE, fp) != NULL) {
        pt = strtok(line, " ");
        
        for (int i = 0; i < sizeof(entries) / sizeof(entries[0]); i++) {
            if (strcmp(pt, entries[i].name) == 0) {
                switch (entries[i].type) {
                    case 'd':
                        pt = strtok(NULL, " ");
                        sscanf(pt, "%d", (int *)entries[i].ptr);
                        break;
                    case 'f':
                        pt = strtok(NULL, " ");
                        sscanf(pt, "%lf", (double *)entries[i].ptr);
                        break;
                    case 's':
                        pt = strtok(NULL, " ");
                        sscanf(pt, "%s", (char *)entries[i].ptr);
                        break;
                    case 'p':
                        printf("Reading %s parameters\n", (char *)entries[i].ptr);
                        break;
                    case 'n':
                        // Just a newline, do nothing
                        break;
                    default:
                        printf("Keyword unknown: %s\n", pt);
                        break;
                }
                break;
            }
        }
    }

    psl->nparts=sys.npart;
    sys.beta=psl->beta;
    fclose(fp);
    
    printf("Done reading path.inp\n");

    return;
}


void init_model(Slice *psl) {
    /*startis with defining the potential based on rwet and sc,
    then determines the minimum of the potential and de bond energy threshold
    */
    double bond_treshold_fraction=0.001;

    sys.sigma=3.2e-6;//micron
    if (sys.boxl.y<1e-5 & sys.boxl.x>1e-5){
        //its a cubic box
        sys.boxl.y = sys.boxl.z = sys.boxl.x;
    }
    else if(sys.boxl.y>1e-5 & sys.boxl.x>1e-5){
        sys.boxl.z = sys.boxl.x;
    }
    else if(sys.boxl.x<1e-5){
        vprint(sys.boxl);
        error("box length is not read correctly");
    }

    psl->temp = 1.0/psl->beta;
    cluster.update=1; // always start with a cluster update if dong MC + clusterupdate 
    pot.rcutoff=pot.s_cutoff+1.;

    /*setting up the parameters of the potentials, based on wetting , surface charge*/
    setup_criticalCasimir_potential_parameters(); // casimir potential 

    /*take minimum bond length for bond average*/
    printf("\nLooking for the minimum of the potential...");
    pot.s_min= find_minimum_of_potential();
    pot.s_overlap = find_overlap_distance(); // do this after finding s_min

    // NOTE: only used in BMD!! MOVE THIS FUNCTION
    // pot.s_forcecut = find_forcecutoff_distance(); 

    pot.bond_cutoffE=-0.1; // this defines the bond based on energy

    //first read in the particle, specifies particle types etc
    printf("Setting up the positions and sites for the particles...");
    setup_positions_sites(psl);
    printf("done\n");
    
    // //patchwidth variables
    // printf("Setting up the patch widths definitions for the particles...");
    // setup_delta();
    // printf("done\n");

    // printf("printing potentials ... \n");
    // plotpotential();
    // printf("done\n");


    // /*Calculate the gravitational force fg and correction b_zc*/
    // if(sys.gravity>0){
    //     printf("\nSetting GRAVITY PARAMETERS...\n");
    //     gravitational_parameters();
    //     printf("done\n");
    //     double cov=psl->nparts/(sys.boxl.x*sys.boxl.y);
    //     printf(" the density of the quasi-2D system %.4lf [N/sigma^2] = %.4lf area percentage ",cov,cov*(PI/4.));
    // }
    
    // /*rcut^2*/
    // if (pot.s_cutoff>0.5){
    //     error("pot.s_cutoff too bit, make it bigger than 0.5");
    // }
    
   
    // /*does sim_type specific definitions*/
    // init_simtype(psl);
    // check_input_with_MAXDEFS();

    return;
}
double find_minimum_of_potential(void){
    /* actually just use deriveative much easier? 
    this code looks for the minimum of the potential by increasing s the surface-surface distance and reading the attractive and repulsive potential
    it returns the surface-surface-distance at the minimum*/

    double s_min = 1.1, Utot_tracked = 0, Uattr_tracked = 0, Urep_tracked = 0;
    double s, Urep, Uattr, Utot;

    // HJ: dont remember why this is necessary.. I think
    potential_attractive_energy_sdist(0.);

    // if you are using the square well potential, return the cutoff (which is an input variable)
    if(pot.sqw_epsilon!=0){ return pot.s_cutoff;}
    else{ // when using critical Casimir
        for (int r=1;r<100000;r++){
            s=r/100000. *0.10;
    
            Urep = potential_repulsive_energy_sdist(s);
            Uattr = potential_attractive_energy_sdist(s);
            Utot = Urep+Uattr;
            
    
            if (Utot<Utot_tracked){
                s_min=s;
                Utot_tracked=Utot;
                Uattr_tracked=Uattr;
                Urep_tracked=Urep;
    
            }
        }
        pot.Erep_smin=Urep_tracked;
        pot.Ec_smin=Uattr_tracked;
        pot.E_smin=Utot_tracked;

        printf("    the colloid-colloid distance with minimum energy sys.E_smin=%lf (sys.Ec_smin=%lf, sys.Erep_smin=%lf) is %lf \n", pot.E_smin, pot.Ec_smin, pot.Erep_smin, s_min );
        pot.s_min=s_min;
        printf("    found s_min = %lf [sigma]\n", pot.s_min);
        
        if (pot.s_min>1.0){
            error("!!define pot.s_min such that the minimum energy and the bond definition ( == 1 procent of minimum ) can be calculated ");
        }
        return s_min;
    }
}
double find_overlap_distance(void ){
    /*looks for the value of s (surface-surface-distance) at which you define the particles to overlap = s_overlap
    take it as the distance at which the  energy is 22/beta away from the minimum  --> probability of 1 in 4 miljard = exp(-22)
    */

    double Delta_E=50./sys.beta,s,Uattr,Urep,Utot,dE ;
    int x=1000000;       
    gprint(pot.E_smin);

    for (int r=x/-5;r<x +1;r++){ //  start at negative values for LJ
        s=(double)r/(double)x *pot.s_cutoff;
        
        Urep = potential_repulsive_energy_sdist(s);
        Uattr = potential_attractive_energy_sdist(s);
        Utot = Urep+Uattr;
        // printf(" at distance s=%3.2f Urep==%3.2f +  Uattr=%3.2f  = %3.2f \n",s,Urep,Uattr,Utot);
        dE=Utot-pot.E_smin; // dE away from the minimum
        if ((Delta_E-dE)>1e-1){
            printf(" at distance s=%3.2f Urep==%3.2f +  Uattr=%3.2f  =  %3.2f\n",s,Urep,Uattr,Utot);

            gprint(Delta_E);
            gprint(dE);
            printf("found the overlap distance  = %lf [sigma]\n", s);    
            return s;

        }
    }

    printf("found the overlap distance  = %lf [sigma]\n", s);
    return s;

}
double find_forcecutoff_distance(void ){
    /*looks for the value of s (surface-surface-distance) at which you the force is large and may causes the particles to explode
    */ 

    double Delta_E=-pot.E_smin/sys.beta,s,Uattr,Urep,Utot,dE ;
    int x=100000;       
    gprint(pot.E_smin);


    for (int r=x/-5;r<x +1;r++){ //  start at negative values for LJ
        s=(double)r/(double)x *pot.s_cutoff;
        
        Urep = potential_repulsive_energy_sdist(s);
        Uattr = potential_attractive_energy_sdist(s);
        Utot = Urep+Uattr;

        dE=Utot-pot.E_smin; // dE away from the minimum
        if (fabs(Delta_E-dE)<1e-1){
            printf(" at distance s=%3.20f Urep==%3.2f +  Uattr=%3.2f  =  %3.2f\n",s,Urep,Uattr,Utot);
            printf(" at distance s=%3.20f Frep==%3.2f +  Fattr=%3.2f  \n",s,bond_repulsive_force(s),bond_attractive_force(s));

            gprint(Delta_E);
            gprint(dE);
            return s;
        }
    }
    // if there is no such distance, just take s_overlap
    return pot.s_overlap;

}

void setup_positions_sites(Slice *psl) {

    double r2,s,radius_i;
    int ipart, jpart, overlap, wall,n;
    vector dr;
    Pts *psi, *psj;
    int maxchain_length, mod_ipart; // Declare the variables here


    read_particletypes(psl);
    // print_particle_properties();

    switch (init_conf.start_type) {
        case 0:   
            printf("    Randomly place the particles with random orientation and position\n");
            for (ipart = 0; ipart < psl->nparts; ipart++) {
                psi = &psl->pts[ipart];
                radius_i = sys.particletype[psi->ptype].radius;
                psi->q = RandomQuaternion();

                do {
                   /*place them randomly; do not put the particle at z=[0,1] if gravity is on*/
                    overlap = 0;
                    psi->r.x = RandomNumberRange(-0.5 * sys.boxl.x, 0.5 * sys.boxl.x);
                    psi->r.y = RandomNumberRange(-0.5 * sys.boxl.y, 0.5 * sys.boxl.y);

                    if (sys.gravity > 0)
                        psi->r.z = RandomNumberRange(1.12 + radius_i, 1.12 + 2. * radius_i);
                    else
                        psi->r.z = RandomNumberRange(-0.5 * sys.boxl.z, 0.5 * sys.boxl.z);
                    
                    // check for overlap
                    for (jpart = 0; jpart < ipart; jpart++) {
                        psj = &psl->pts[jpart];
                        vector_minus(psj->r, psi->r, dr);
                        pbc(dr, sys.boxl);
                        r2 = vector_inp(dr, dr);

                        s = sqrt(r2) - radius_i - sys.particletype[psj->ptype].radius;
                        if (s < pot.s_min * 0.9)
                            overlap = 1;
                    }
                } while (overlap == 1);
            }
            break;
        case 1:   
            /*read from input*/
            if(init_conf.start_type==1) {
                    printf("Reading from conf.inp\n");
                    conf_input(&slice[0]);
            }
            break;
        case 2:
            /* place them in a chain*/
            
            printf("Placing %d particles in single chain\n", psl->nparts);
            double y_previousparticle=-0.48*sys.boxl.y;
            if((double)sys.boxl.x/psl->nparts<=pot.s_min){
                error("ERROR: Too many particles to put in a 1D chain. Adjust boxl or npart!");
            }
            else{
                printf("pot.s_min = %lf\n", pot.s_min);      
                for( ipart=0; ipart<psl->nparts; ipart++) {
                    psi=&psl->pts[ipart]; 
                    psi->r=nulvec;
                    psi->r.z =1.15+sys.particletype[psi->ptype].radius;
                    psi->r.y = y_previousparticle + (sys.particletype[psi->ptype].radius+pot.s_min) ;                    
                    psi->q.q0  =1.;
                    psi->q.q1 = psi->q.q2 = psi->q.q3 = 0.;
                    y_previousparticle = psi->r.y+sys.particletype[psi->ptype].radius;
                }
            }  
            break;
        case 4:
            // first check if your chains will fit the box
            maxchain_length=floor(psl->nparts/init_conf.nchains)+1;

            printf("Placing %d particles in %d chains\n", psl->nparts,init_conf.nchains);
            if((double)(sys.boxl.y<=maxchain_length*(pot.s_min+1) )|| (double)(sys.boxl.x)<=init_conf.chaingap*(init_conf.nchains)){
                printf("ERROR: Too many particles to put in a 1D chain %lf<%lf or %lf<%lf . Adjust boxl or npart! \n", (double)sys.boxl.y,maxchain_length*(pot.s_min+1),sys.boxl.x,init_conf.chaingap*(init_conf.nchains+1));
                exit(1);
            }
            else{
                n=0;

                for( ipart=0; ipart<psl->nparts; ipart++) {
                    mod_ipart=ipart%(maxchain_length+1);
                    psi=&psl->pts[ipart]; 

                    psi->r.x = -0.48*sys.boxl.x + init_conf.chaingap*n;
                    psi->r.z = 1.13+sys.particletype[psi->ptype].radius;
                    psi->r.y = -0.48*sys.boxl.y + (pot.s_min+2.*(sys.particletype[psi->ptype].radius))*mod_ipart ;   
                    
                    psi->q.q0=1.;
                    psi->q.q1=0.;
                    psi->q.q2=0.;
                    psi->q.q3=0.;
                    
                    if(ipart%maxchain_length==0 & ipart>0) n++;
                    
                }  
            }  
            break;
    }
    
    
    printf("Setting the positions for all the particles is done\n");

    /*setup al the patch vectors*/
    update_patch_vectors(psl);

    
    return;
}

void conf_input(Slice *psl) {

    int ipart, isite,npart,nsites;
    vector boxl;
    FILE *fp;

    if((fp=fopen("conf.inp","r"))==NULL) {;
        printf("Warning: can not open conf.inp\n");
    }
    else {
        fscanf(fp,"%d %lf %lf %lf\n", &npart, &boxl.x, &boxl.y, &boxl.z);
        for(ipart=0; ipart<npart; ipart++) {
            fscanf(fp,"%lf %lf %lf ", &psl->pts[ipart].r.x, &psl->pts[ipart].r.y, &psl->pts[ipart].r.z);
            fscanf(fp,"%lf %lf %lf %lf\n", &psl->pts[ipart].q.q0, &psl->pts[ipart].q.q1, &psl->pts[ipart].q.q2, &psl->pts[ipart].q.q3);

             /*make sure the quaternions are unit, it could be a liiiiilte bit off due to the print accuracy and it will  lead to |p|!=1 */
            psl->pts[ipart].q=  normalize_quat(psl->pts[ipart].q);
            pbc(psl->pts[ipart].r,sys.boxl);
        }
    }
    fclose(fp);

    if(npart!=psl->nparts) {
        error("Error: number of particles in system not same as in conf.inp\n");
        
    }
    if((boxl.x!=sys.boxl.x) | (boxl.y!=sys.boxl.y) | (boxl.z!=sys.boxl.z)) {
        vprint(boxl);
        vprint(sys.boxl);
        error("Error: boxlength in system not same as in conf.inp\n");
    }

    return;
}

// double find_trunc_of_Saccent_1( int ptype){
//     /* find the angle at which S' is almost zero,i.e. < treshold in an iterative manner*/
//     double theta, costheta,S_new,treshold=1e-6,max_angle=30.0;
//     int i,imax=1000000;

//     for(i=0;i<imax; i++){
//         theta=(double)i/imax*max_angle;
//         costheta=cos(theta/180.*PI);        
//         S_new=Saccent(costheta,ptype);

//         if (S_new<treshold){
//             return costheta;
//         }
//     }
    
//     error("cosdelta not found.");
    
//     return -1;
// }

// double find_zcut(double fg){
//     /*fg_LJ = 4.*pot.epsilongravLJ*(-12.*zcutinv12*zcutinv+6.*zcutinv6*zcutinv);
//     zcut between 1 and 1.2*/

//     int z;
//     double z_search,zmax=1e8;
//     double fg_LJ;
//     double zcutinv, zcutinv2,zcutinv6,zcutinv12;

//     for(z=0;z<(int)zmax;z++){
//         z_search = (double)z/zmax*0.2+1.;
//         zcutinv = 1./z_search;
//         zcutinv2 = zcutinv*zcutinv;
//         zcutinv6 = zcutinv2*zcutinv2*zcutinv2;
//         zcutinv12 = zcutinv6*zcutinv6;
//         fg_LJ = 4.*pot.epsilongravLJ*(-12.*zcutinv12*zcutinv+6.*zcutinv6*zcutinv);
//         if(fabs(fg - fg_LJ)<=1e-4){
//             gprint(fg);
//             gprint(fg_LJ);
//             gprint(z_search);
//             return z_search;
//         }
//     }
//     error("zcut not found");
//     return 0;
// }

// void gravitational_parameters(void){
//     double g,delta_rho_kg_m3,fg_LJ, boltzmann, T, ktsigma;
//     double zcutinv, zcutinv2,zcutinv6,zcutinv12;
//     Particletype *ptypen;
//     double fg,zcut,b_zc,r_micron;
//     int n;
    
    
//     // constants
//     g=9.80665; /*gravitation acceleration of eath [m/s^2]*/
//     boltzmann=1.38064852e-5; /*joule. the e-5 comes from e-23 *(original value of Boltzmann) times (e-6)^3=e-18 (because in sys.fg we use sys.r_micron^3) */
//     T=33.8+273.15;/*K, this an is approximate for all measurements*/
//     // already defined sys.sigma=3.2e-6;/*sigma is reducing measure for distance 3.2 mum= the colloidal diameter*/

//     //some particle indep paramters
//     ktsigma = (boltzmann * T)/sys.sigma;
//     printf("boltzmann %.12lf [e-18 joule], T %.12lf [K] ,sigma %.12lf [m]\n", boltzmann,T,sys.sigma);
//     printf("kt/sigma %.12lf \n",ktsigma);
    
//     //some particle dep paramters
//     for (n=0;n<sys.nparticle_types;n++){
//         ptypen=&sys.particletype[n];
//         if((ptypen->nsites==2 )|| (ptypen->nsites==1)){
//             ptypen->delta_rho_kg_m3=60.568460120322576 ; /* the density diffence between the lutide-water solution and the colloid [kg/m^3]*/
//         } 
//         else if((ptypen->nsites==4 )|| (ptypen->nsites==3)){
//             ptypen->delta_rho_kg_m3=63.5699387858425;
//         }
//         else if(ptypen->nsites==NSITES){ //isotropic particle
//             double rho_mixture,rho_TPM,delta_rho;
//             rho_mixture=0.98965956933887;//[g/mL] J. Chem. Phys., Vol. 100, No.1, 1 January 1994 table III
//             rho_TPM=1.235; // patch material [g/mL] from (non polymerized, see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5569361/)
//             ptypen->delta_rho_kg_m3= (rho_TPM-rho_mixture)*1000.; //[g/mL]*1000 = [kg/m^3]
//         }
//         else{
//             dprint(ptypen->nsites);
//             error("npatch not corresponding to 2 or 4. Indicate delta_rho_kg_m3 in the code");
//         }
//     }

//     //calculate the gravity per particle, 
//     for (n=0;n<sys.nparticle_types;n++){
//         ptypen=&sys.particletype[n];
//         //the radius of the particle in micron; sigma is given in meter
//         r_micron=ptypen->radius*sys.sigma*1e6;
        

//         fg = 4./3.*PI* (r_micron*r_micron*r_micron)*ptypen->delta_rho_kg_m3*g/(ktsigma)*sys.gravity; /* V_g = fg*z [kT/sigma*] at 307K*/
//         printf("sys.sigma=%.5lf ,r_micron= %.5lf, fg = %lf [kT/sigma]\n", sys.sigma,r_micron,fg);
//         zcut=find_zcut(fg);
//         zcutinv = 1./(zcut);
//         zcutinv2 = zcutinv*zcutinv;
//         zcutinv6 = zcutinv2*zcutinv2*zcutinv2;
//         zcutinv12 = zcutinv6*zcutinv6;

//         fg_LJ = 4.*pot.epsilongravLJ*(-12.*zcutinv12*zcutinv+6.*zcutinv6*zcutinv);
//         if(fabs(fg - fg_LJ)>=0.01){
//             printf("the difference of Fg(line)-fgLJ = %lf -%lf  = %lf", fg , fg_LJ, fg - fg_LJ);
//             error("gravity force not equal at zcut");
//         }
//         /* b_zc should include correction of radius. The wall is at 0 kT and the center of the colloid is then at z=radius+1.12 (1.12 is minimum of LJ) */
//         b_zc = fg*(zcut) - 4.*pot.epsilongravLJ*(zcutinv12-zcutinv6+1./4) ;
        
//         ptypen->fg=fg;
//         ptypen->zcut=zcut;
//         ptypen->b_zc=b_zc;

//     }

//     return;
// }





// void setup_delta(void){
//     // walk over the particle types and identify what kind of switch function is has, and what the delta cutoff is (i.e. where is S' 0)
//     int ptype;
//     Particletype *part_type;
//     Site *site_type;
//     double mini_costheta=1.;

//     for(ptype=0;ptype<sys.nparticle_types; ptype++){
//         part_type=&sys.particletype[ptype];
//         site_type=&site[ptype];

//      // if delta if >0 then use old switch as S'
//         if (site_type->s_accent==0 ){
//             printf("OLD SWITCH\n"); // (smooth)
//             site_type->cosdelta=cos(site_type->delta_degree/180.*PI); // theta cutoff
//             site_type->oneover_cosdelta = 1.0/(1.0-site_type->cosdelta);
//         }
//         else if ( site_type->s_accent==2 ){
//             printf("LINEAR SWITCH\n");
//             site_type->cosdelta=cos(site_type->delta_degree/180.*PI); // theta cutoff
//         }
//         //if delta<0 (<1e-5) then use integrated switch function
//         else if (site_type->s_accent==1){ // ctiyicl casimir
//             printf("NEW SWITCH\n");
//             site_type->cosdelta=find_trunc_of_Saccent_1(ptype);
//             site_type->delta_degree= acos(site_type->cosdelta)/PI*180.;
//         }
//         else if(site_type->s_accent==3 && site_type->S_fixed>0){
//             site_type->cosdelta=cos(site_type->delta_degree/180.*PI);
//         }
//         else{
//             error("something went wrong in setting up delta (the tresholds of the patch width in e.g. energy calculation.");
//         }

//         if (site_type->cosdelta<mini_costheta){
//             mini_costheta = site_type->cosdelta;
//         }
//     }

//     return;
// }

// void setup_mc_move(MC *mctype){
//         printf("*****setting up the MC %s ****\n",mctype->name);
//         mctype->rot.acc=0;
//         mctype->trans.acc=0;
//         mctype->fintrans.acc=0;
//         mctype->finrot.acc=0;
//         mctype->trans.tries=0;
//         mctype->fintrans.tries=0;
//         mctype->rot.tries=0;
//         mctype->finrot.tries=0;

//         mctype->drmax=1.;
//         mctype->dqmax=5.;
//         return;
// }
// void init_simtype(Slice *psl){
//     int  id;
//     /*if bmd or ffs*/

//     // when restarting the simulation, you also need to read in a configuration from file
//     if ( (sys.restart==1) && (init_conf.start_type!=1)){ 
//         error("you are usign a restart time, but init_conf.start_type!=1. contradicting settings"); 
//     }
//     if (sys.sim_type==SIM_BMD || sys.sim_type==4){
//         //do here the force check, you need s_min
//         printf("\nChecking the derivatives/forces...");
//         derivative_check(psl);
//         printf("done\n");
//         /*langevin/brownian md parameters*/
//         printf("\nDefining the langevin timestep and mobility parameters... ");

//         // so the input used to read mobility and beta, but its odd. 
//         //I changed it to reading in diffusion 13-okt, so mobility is calculated inside the code via beta
//         //https://en.wikipedia.org/wiki/Einstein_relation_(kinetic_theory)
//         // D = mu*kT = mu/Beta ; thus mu = D*beta
//         if (langevin.mobilityT <1e-5 || langevin.mobilityR <1e-5){
//             gprint(langevin.mobilityT);
//             gprint(langevin.mobilityR);
//             error(" the rotational diffusion and/or translational diffusion/mobiilty is close to zero? ");
//         }
//         langevin.diffusionT= langevin.mobilityT/sys.beta;
//         langevin.diffusionR= langevin.mobilityR/sys.beta;
        
//         langevin.sqrtmobilityT = sqrt(langevin.mobilityT); //eqn 5,  
//         langevin.sqrtmobilityR = sqrt(langevin.mobilityR); // eqn 4 Ilie2015

//         langevin.dtD=sqrt(2.0*langevin.timestep*psl->temp);
//         //langevin.dtBeta=langevin.timestep*psl->beta; // oeeehhh I see that beta is put in here... 

        
//         // depending on the input you give, you determine the ncycle1, langevin.step
//         if (sys.ncycle1==0){
//             langevin.step=(int)ceil(langevin.print_time/(langevin.timestep*langevin.ninter));
//             sys.ncycle1=langevin.total_time/(langevin.step*langevin.timestep*langevin.ninter*sys.ncycle2);
//         }
//         else{
//             langevin.step=100;
//         }

//         if(sys.restart){
//             restart_read_time(psl);
//         }

//         if (sys.ncycle1==0 || sys.ncycle2==0 || langevin.step==0 || langevin.ninter==0){
//             dprint(sys.ncycle1);
//             dprint(sys.ncycle2);
//             dprint(langevin.step);
//             dprint(langevin.ninter);
//             error(" one of these is zero. stop");
//         }
        
//         if (sys.mc_warmup>0){
//             sprintf(mc_single_large.name,"single particle large moves");
//             setup_mc_move(&mc_single_large);

//             sprintf(mc_single_small.name,"single particle small moves");
//             setup_mc_move(&mc_single_small);

//         }

//         printf("done\n");

//             /* create the neighbor list here;*/
//         if(sys.nearest_neighbor==1){
//             printf("Setting up the neighbor list\n");
//             setup_nnlist();
//             printf("making the neighbor\n");
//             update_nnlist(psl);
//         }
//     }
//     else if(sys.sim_type==SIM_MC) {
//         //Monte Carlo
//         // printf("\nSetting up the MC parameters... ");
//         //Monte Carlo; there is only MC in this code...

//         if (sys.nearest_neighbor==1 & sys.sim_type==SIM_MC ){
//             error("nearest_neighbor and MC do not go together!");
//         }

//         psl_old=(Slice *)calloc(1,sizeof(Slice)); //global variable,used in cluster MC
//         copyslice=(Slice *)calloc(1,sizeof(Slice)); //global variable,used in cluster MC
//         printf("\nSetting up the MC parameters...\n ");

//         sprintf(mc_single_large.name,"single particle large moves");
//         setup_mc_move(&mc_single_large);

//         sprintf(mc_single_small.name,"single particle small moves");
//         setup_mc_move(&mc_single_small);
//         mc_single_small.drmax/=100.;
//         mc_single_small.dqmax/=10.;

       

//         gprint(mc_single_large.drmax);
//         gprint(mc_single_large.dqmax);
        
//         if(sys.cluster_MC==1){
//             sprintf(mc_cluster.name,"cluster moves");
//             setup_mc_move(&mc_cluster);

//             sprintf(mc_mono.name,"single-particle-clusters moves");
//             setup_mc_move(&mc_mono);
//         }  
//         printf("   .. done\n");  
//     }


//     if(cluster.analysis==1){
//         // chain. = (Statistics *)calloc(NPART,sizeof(Statistics));
//         printf("\nCode performs cluster analysis \n");
//         printf("the total energy is %lf\n", total_energy(&slice[0]));
//         printf("now cluster analysis\n");
//         slice[0].nclusters = cluster_analysis(&slice[0]); 
//         clustersize_identification(&slice[0]);

//         printf("Initial # bonds is %d \n", slice[0].nbonds);
     
//         strcpy(cluster.size_histogram.filename,"clustersize_histogram.out");
//         strcpy(cluster.size_distribution.filename,"clustersize_distribution.out");

//         if(sys.restart){
//             read_statistics_file(&cluster.size_histogram);
//             read_statistics_file(&cluster.size_distribution);          
//         }   
//     }

//     return;
// }
// void restart_read_time(Slice *psl){
//     // read the last line of the trajectory.xyz (just copy last line to a new file)
//     // code from https://stackoverflow.com/questions/13790662/c-read-only-last-line-of-a-file-no-loops

//     FILE *fd;                               // File pointer
//     char *pt;                         
//     char filename[] = "trajectory.xyz";       // file to read
//     static const long max_len = 150+ 1;  // define the max length of the line to read
//     char buff[max_len + 1];             // define the buffer and allocate the length
//     double dummy;
//     // printf(">>restart_read_time<<\n");
//     if ((fd = fopen(filename, "rb")) != NULL)  {      // open file. I omit error checks

//         fseek(fd, -max_len, SEEK_END);            // set pointer to the end of file minus the length you need. Presumably there can be more than one new line caracter
//         fread(buff, max_len-1, 1, fd);            // read the contents of the file starting from where fseek() positioned us
//         fclose(fd);                               // close the file

//         buff[max_len-1] = '\0';                   // close the string
//         char *last_newline = strrchr(buff, '\n'); // find last occurrence of newlinw 
//         char *last_line = last_newline +1;         // jump to it

//         // printf("captured: [%s]\n", last_line);    // captured: [472977827]

//         pt = strtok(last_line," ");
//         sscanf(pt,"%lf",&psl->c_time);
//         printf("  restart time = %lf [s]\n", psl->c_time);    // captured: [472977827]
//     }
//     else{
//         psl->c_time=0.;
//         printf(" >>>>>>>WARNING:  trajectory.xyz cannot be read , so restart time = %lf [s]<<<<<\n", psl->c_time); 

//     }

// }
// void plotpotential(void){
//     FILE *fp;
//     int i,j,n, t;
//     double r,s,z, V=0 , V_LJ, rinv, rinv3, rinv6,rinv12,rinv24, V_g, Vrep=0, Vattr=0;
//     double Eattr=0, Erep=0, E=0;
//     double theta, costheta,S_old, S_new,S,dS,cosdelta_n; 
    
//     int x=1000;

//     if ((fp = fopen("pot.dat","w"))==NULL){
//         printf("output: can't open pot.dat\n");
//         return;
//     }
//     else {
   
//         for(i=x/-5;i<x+1; i++){
//             s=(double)i/(double)x*pot.s_cutoff;
            
//             Vrep = potential_repulsive_energy_sdist(s);
//             Vattr = potential_attractive_energy_sdist(s);
            
//             V= (Vrep + Vattr);
            
//             fprintf(fp, "%30.15lf %30.15lf %30.15lf %30.15lf\n", s, V, Vrep, Vattr);//, E, Erep, Eattr);
//         }
//     fclose(fp);
//     }
//     if (sys.sim_type==0){ 
//         if ((fp = fopen("2nd_der.dat","w"))==NULL){
//             printf("output: can't open 2nd_der.dat\n");
//             return;
//         }
//         else {
//             // printf("    note: 2nd_der.dat does not include Derjaguin scaling\n");

//             for(i=x/-5;i<x+1; i++){
//                 s=(double)i/(double)x*pot.s_cutoff;
                
//                 Vrep = second_der_Vrep(s);
//                 Vattr = second_der_Vc(s);
                          
//                 V= Vrep + Vattr;
                
//                 fprintf(fp, "%30.15lf %30.15lf %30.15lf %30.15lf\n", s, V, Vrep, Vattr);//, E, Erep, Eattr);
//             }
//         fclose(fp);
//         }
    
//         if ((fp = fopen("force.dat","w"))==NULL){
//             printf("output: can't open force.dat\n");
//             return;
//         }
//         else {
//             // printf("    note: force.dat does not include Derjaguin scaling\n");

//             for(i=x/-5;i<x+1; i++){
//                 s=(double)i/(double)x*pot.s_cutoff;
                
//                 Vrep = bond_repulsive_force(s);
//                 Vattr = bond_attractive_force(s);
                     
//                 V= Vrep + Vattr;
                
//                 fprintf(fp, "%30.15lf %30.15lf %30.15lf %30.15lf\n", s, V, Vrep, Vattr);//, E, Erep, Eattr);
//             }
//         fclose(fp);
//         }
//     }


//     if ((fp = fopen("switchfunction.dat","w"))==NULL){
//         printf("output: can't open switchfunction.dat\n");
//         return;
//     }
//     else {
        
//         for(i=0;i<1000; i++){
//             theta=i/1000.*50.0;
//             costheta=cos(theta/180.*PI);
//             fprintf(fp, "%30.15lf ",theta);
//             for(n=0;n<sys.nparticle_types;n++){
//                 cosdelta_n=site[n].cosdelta;

//                 if(costheta<cosdelta_n) {
//                     S=0;
//                     dS=0;
//                 }
//                 else{
//                     S=Saccent(costheta,n); 
//                     if (sys.sim_type!=2){
//                         dS=dSaccent_dcosangle(costheta, n);
//                     }
//                     else{
//                         dS=0;
//                     }
//                 }   
//                 fprintf(fp, "%30.15lf %30.15lf",S,dS);
//             }
            
//             fprintf(fp, "\n");//, E, Erep, Eattr);
//         }
//     fclose(fp);
//     }


//     if(sys.gravity>0){
//         if ((fp = fopen("gravity.dat","w"))==NULL){
//             printf("output: can't open gravity.dat\n");
//             return;
//         }
//         else {
            
//             for(i=100;i<1000; i++){
//                 r=i/1000.*3+1;
                
//                 fprintf(fp, "%30.15lf ", r);

//                 for(n=0;n<sys.nparticle_types;n++){
//                     z=r-sys.particletype[n].radius;

//                     if(z>=sys.particletype[n].zcut){ 
//                         V_g = sys.particletype[n].fg*z - sys.particletype[n].b_zc ; /*V_g = Fg*z - b*/
//                     }
//                     else{
//                         rinv = 1/z;
//                         rinv3= rinv*rinv*rinv;
//                         rinv6=rinv3*rinv3;
//                         rinv12= rinv6* rinv6;
//                         rinv24 = rinv12*rinv12;
//                         V_g = 4*pot.epsilongravLJ*(rinv12-rinv6+1./4); /* LJ12-6 = 4*b(1/r12-1/r6)*/
//                     }

//                     fprintf(fp, "%30.15lf ", V_g);
//                 }
//                 fprintf(fp, "\n");
//             }
            
//         fclose(fp);
//         }
//     }


//     return;
// }
   

// void check_input_with_MAXDEFS(){
//     /*in path.h there are max definitions given. 
//     Check after reading the input file if any of these max values are crossed*/

//     if (sys.npart>NPART){
//         dprint(sys.npart);
//         dprint(NPART);
//         error("sys.npart > NPART");
//     }
//     if (sys.nparticle_types>PTYPES){
//         dprint(sys.nparticle_types);
//         dprint(PTYPES);
//         error("sys.nparticle_types > PTYPES");
//     }
//     for (int p=0; p<sys.nparticle_types;p++){
//         if(sys.particletype[p].nsites>NSITES){
//             error("part_type->nsites > NSITES");
//         }
//     }
    

//     return;

// }



void read_saccent_sfixed_from_pathinp(void){
    /*the s_fixed, s_accent are read from path.inp*/
    char *pt,line[MAXLINE];
    FILE *fp;
    char filename[100];
    Site *st;
    Particletype *pt_type;

    if((fp = fopen("path.inp","r"))==NULL) {
        printf("ERROR: could not read input file\n");
        exit(1);
    }
    
    printf("  reading from path.inp: \n");
   
    while(fgets(line,MAXLINE, fp) != NULL) {
        pt = strtok(line," ");
        // printf("%s\n",line);
        if( strcmp(pt,"s_accent")==0) {  
           printf("      s_accent  "); 
           int n=0;
           pt = strtok(NULL," ");
           // walk through other tokens 
           while( pt != NULL ) {
                st=&site[n];
                sscanf(pt,"%d",&st->s_accent);  // read all s_accents
                pt = strtok(NULL, " ");
                n++;
                if( st->s_accent>3 || st->s_accent<0){
                    dprint(st->s_accent);
                    error("st->s_accent can only  be 0, 1, 2, or 3 (resp. old, integrated, linear, fixed)  ");
                }
           }
           if (n!=sys.nparticle_types){
                dprint(n);
                dprint(sys.nparticle_types);
                printf("mind that there cannot be a space behind the variables, that may cause this error: \n");
                error("too many s_accents in path.inp");
           }
           printf(" ... done \n");
        }
        else if( strcmp(pt,"S_fixed")==0) {   
           printf("      S_fixed "); 
           int n=0;
           pt = strtok(NULL," ");
           // walk through other tokens 
           while( pt != NULL ) {
                st=&site[n];
            
                sscanf(pt,"%lf",&st->S_fixed);  // read all S_fixeds 
                pt = strtok(NULL, " ");
                // printf("reading from path.inp the S_fixed n=%d\n",n);

                n++;
                if( st->S_fixed>1 || st->S_fixed<0){
                    gprint(st->S_fixed);
                    error("st->S_fixed can only be between [0,1] ");
                }
           }
           if (n!=sys.nparticle_types){
                dprint(n);
                dprint(sys.nparticle_types);
                printf("mind that there cannot be a space behind the variables, that may cause this error: \n");
                error("too many S_fixeds in path.inp");
           }
           printf(" ... done \n");

        }
    } 
     
    fclose(fp);
    return;

}
void read_particletypes(Slice *psl){
    /*read sites and diameter from particle.inp (used to be sites.inp)*/
    /* the file is structured as follows:
    NPARTICLE_TYPES
    NSITES NPARTICLES DIAMETER DELTA ACTIVITY
    SITE_1 
    SITE_N
    ACTIVITY_VECTOR
    etc.
    */
    FILE *file;
    int isite, ptype,a,p,ipart=0;
    int amount, totalparticles=0;
    vector new;
    Particletype *part_type;
    Site *st;
    double activity_p;
    

    
    if(init_conf.particle_setup==1) {
        printf("Reading particles setup directly from particles.inp\n");
        
        if ((file = fopen("particles.inp","r"))==NULL){
            error("input: can't open particles.inp \n");
        }
        else {
            /*first line of sites.inp contains the number of particle types*/
            fscanf(file,"%d\n",&sys.nparticle_types);
            if (sys.nparticle_types>PTYPES){
                dprint(sys.nparticle_types);
                dprint(PTYPES);
                error("there is a boundary of maximum number of particle types");
            }

            /*such that you know how many to read in*/
            for(ptype=0;ptype<sys.nparticle_types; ptype++){
                /*the first line of the particle type contains:
                number of sites, number of particles, diameter in sigma*/
                part_type=&sys.particletype[ptype];
                st=&site[ptype];

                // not put s_accent, and Rp in particles.inp because then the code becomes different form the old format. which might lead to accidental bugs
                // fscanf(file,"\n%d %d %d %lf %lf %lf\n",&part_type->nsites, &part_type->nparticles, &st->s_accent, &part_type->diameter, &st->delta_degree, &activity_p);
                fscanf(file,"\n%d %d %lf %lf %lf\n",&part_type->nsites, &part_type->nparticles, &part_type->diameter, &st->delta_degree, &activity_p);
                part_type->radius=part_type->diameter/2.;


                /*track if all particles have an assigned patch/diameter definition*/
                amount =  part_type->nparticles;
                totalparticles += part_type->nparticles;
                sys.npatch+=(amount*part_type->nsites);

                /*read the sites, and make sure these vectors are unit vectors*/
                for(isite=0; isite<part_type->nsites; isite++) {
                    fscanf(file,"%lf %lf %lf\n",
                            &new.x,
                            &new.y,
                            &new.z);
                    part_type->site[isite]= check_read_unitvec(new);
                }
                
                if (activity_p>0){
                    part_type->activity=1;
                    part_type->F_A=activity_p;
                    fscanf(file,"%lf %lf %lf\n",
                            &new.x,
                            &new.y,
                            &new.z);
                    part_type->e_A=check_read_unitvec(new); 
                }
                do{ /*give the individual particles in the slice a particle type*/
                    psl->pts[ipart].ptype = ptype;

                    amount--;
                    ipart++;
                }
                while(amount>0 && ipart<psl->nparts);
            }
            /* check if all particles have an assigned structure (patches and diameter)*/
            if(totalparticles!=sys.npart){
                dprint(totalparticles);
                dprint(psl->nparts);
                dprint(sys.npart);
                error("totalparticles  !=psl->nparts  ");
            }
            fclose(file);
        }

        /* (31-jan_22) in path.inp there is the s_accent listed; its there because of history of the code.
         in the previous version of particles.inp, there was no info on s_accent, so I kept it like that to avoid bugs
        I recommend to change the code if somebody new is going to use the code. Put S_accent in particles.inp */
        read_saccent_sfixed_from_pathinp();

        // read the characteristics of the site, e.g. coefficients of switch function; only when s_accent==1
        for(ptype=0;ptype<sys.nparticle_types; ptype++){

            st=&site[ptype];
            if(st->s_accent!=1){
                continue;
            }
            char filename[100];
            snprintf(   filename, sizeof( filename ), "Scoefficients%d.inp", ptype );
            if ((file = fopen(filename,"r"))==NULL){
                printf("input: can't open %s \n",filename);
                error("input: can't open the file \n");
            }
            else {
                //first indicate if you have Sfixed, Sold or Snew 
                char *pt,line[MAXLINE];
                while(fgets(line,MAXLINE, file) != NULL) {
                    pt = strtok(line," ");
                    if( strcmp(pt,"switch_expc")==0) {   
                        pt = strtok(NULL," ");
                        sscanf(pt,"%lf",&st->s_int.switch_expc);
                    } else if( strcmp(pt,"switch_expd")==0) {   
                        pt = strtok(NULL," ");
                        sscanf(pt,"%lf",&st->s_int.switch_expd);
                    } else if( strcmp(pt,"switch_expe")==0) {   
                        pt = strtok(NULL," ");
                        sscanf(pt,"%lf",&st->s_int.switch_expe);
                    } else if( strcmp(pt,"switch_expf")==0) {   
                        pt = strtok(NULL," ");
                        sscanf(pt,"%lf",&st->s_int.switch_expf);
                    } else if( strcmp(pt,"switch_expg")==0) {   
                        pt = strtok(NULL," ");
                        sscanf(pt,"%lf",&st->s_int.switch_expg);
                    } else if( strcmp(pt,"switch_exph")==0) {   
                        pt = strtok(NULL," ");
                        sscanf(pt,"%lf",&st->s_int.switch_exph);
                    } else if( strcmp(pt,"switch_expi")==0) {   
                        pt = strtok(NULL," ");
                        sscanf(pt,"%lf",&st->s_int.switch_expi);
                    } 
                    else if( strcmp(pt,"switch_unit")==0) {   
                        pt = strtok(NULL," ");
                        sscanf(pt,"%d",&st->s_int.switch_unit);
                    } 
                }               
                fclose(file);
            }
        }
        
    }
    else{
        error("init_conf.particle_setup can only be 1");
    }

    return;
}

vector check_read_unitvec(vector source){
    /*used when reading in vectors from particle.inp, check if the vector is a finite size (>0).
    it normalizes the vector before saving it*/
    /*if the source vector is super small e.g. (0,0,0), it is probably an error*/
    vector source_new;
    if ( fabs(source.x) <1e-5 && fabs(source.y) <1e-5 && fabs(source.z)<1e-5){
        gprint(source.x);
        gprint(source.y);
        gprint(source.z);
        error("patch site not well defined.");
    }
    /*normalize and return*/
    normvec(source,source_new);

    return source_new;
}

// void read_statistics_file(StatsLength *stats_name){
//     FILE *file;
//     char *pt,line[NPART], filename[100];
//     int i=0,j=0, dummy;

//     memcpy(filename,stats_name->filename,sizeof(filename));
//     // printf("filname %s\n", filename);
//     if ((file = fopen(filename,"r"))==NULL){
//         printf("%s can't be opened,  statistics are zero\n", filename);
//     }
//     else{
//         printf("reading %s ...", filename);

//         while(fgets(line,1000, file) != NULL) {
//             sscanf(line,"%d %lf %lf %ld", &dummy, &stats_name->bin[i].mean, &stats_name->bin[i].variance2,&stats_name->bin[i].n);
//             i++;
//         }
//         printf("done\n");
//     }
   
//     return;
// }




// void check_random(void){
//     double number;
//     int bin,nbins=20, histogram[20]={0};
//     double dbin=1./((double )nbins);
//     gprint(dbin);


//     for (int i=0 ; i<1e6 ; i++){
//         number=RandomNumber();
//         bin=(int)(number/dbin);
//         histogram[bin]+=1;
//         // if (i==0){
//         //     gprint(number);
//         //     dprint(bin);
//         // }
//     }

//     FILE *fp;
    
//     if ((fp = fopen("random_hist.dat","w"))==NULL){
//         printf("output: can't open random_hist.dat\n");
//         return;
//     }
//     else {
        
//         for (int i=0 ; i<nbins ; i++){
//             fprintf(fp, "%d %lf\n", histogram[i], (double)histogram[i]/1e6);//, E, Erep, Eattr);
//         }
//     fclose(fp);
//     }
//     return;


// }
// void print_input(Slice *psl) {
    
    
//     printf("Simulation\n");
//     if (sys.read_path){
//         printf("********* READING FROM TRAJECTORYFILE *********\n");
//         printf("directorypath       %s\n", sys.directorypath);


//     }
//     if(sys.sim_type==0){
//         printf("********* BMD *********\n");
//         printf("mobilityT       %.15lf\n", langevin.mobilityT);
//         printf("mobilityR       %.15lf\n", langevin.mobilityR);
//         printf("timestep        %.1lf [ns] / %.1lf [ps]\n", langevin.timestep*1e9, langevin.timestep*1e12);
//         printf("ninter           %d\n", langevin.ninter);
//         printf("nearest_neighbor %d\n", sys.nearest_neighbor);
//     }
//     else if(sys.sim_type==1){
//         printf("********* TIS + BD *********\n");
//         printf("mobilityT       %.15lf\n", langevin.mobilityT);
//         printf("mobilityR       %.15lf\n", langevin.mobilityR);
//         printf("timestep        %.1lf [ns] / %.1lf [ps]\n", langevin.timestep*1e9, langevin.timestep*1e12);
//         printf("ninter           %d\n", langevin.ninter);
//         printf("nearest_neighbor %d\n", sys.nearest_neighbor);
//     }
//     else if(sys.sim_type==2){
//         printf("********* MC *********\n");
//         printf("cluster_MC %d\n", sys.cluster_MC);
        

//     }
//     // else if(sys.sim_type==4){
//     //     printf("********* FFS *********\n");
//     //     printf("mobilityT       %.15lf\n", sys.mobilityT);
//     //     printf("mobilityR       %.15lf\n", sys.mobilityR);
//     //     printf("timestep        %.1lf [ns] / %.1lf [ps]\n", langevin.timestep*1e9, langevin.timestep*1e12);
//     //     printf("ninter           %d\n", langevin.ninter);
//     //     printf("nearest_neighbor %d\n", sys.nearest_neighbor);
//     // }


//     printf("ncycle1          %d\n", sys.ncycle1);
//     printf("ncycle2          %d\n", sys.ncycle2);


//     printf("\nSystem\n");
//     printf("mc_warmup         %d \n", sys.mc_warmup);
//     printf("npart             %d\n", psl->nparts);
//     printf("#particle types   %d\n", sys.nparticle_types);
//     printf("beta              %lf\n", psl->beta);
//     printf("boxl.x            %lf\n", sys.boxl.x);
//     printf("boxl.y            %lf\n", sys.boxl.y);
//     printf("boxl.z            %lf\n", sys.boxl.z);

//     if(init_conf.start_type==2 ||init_conf.start_type==4 ){
//         printf("nchains           %d\n", init_conf.nchains);
//         printf("chaingap x        %lf\n", init_conf.chaingap);
//     }

//     printf("\nPotential Simon\n");
//     printf("rwet              %.3lf\n", pot.r_wetting);
//     printf("surface_charge    %.3lf\n", pot.surface_charge);
//     printf("dl                %.12lf\n", pot.dl);
//     printf("dT                %lf\n", pot.dT);
//     printf("xi(dT)            %.18lf\n", pot.xi);
//     printf("A(dT)             %.18lf\n", pot.A);
//     printf("Arep              %.18lf\n", pot.Arep);
//     printf("s_cutoff           %.3lf\n", pot.s_cutoff);
//     printf("gravity           %lf\n", sys.gravity);
//     if (pot.wall_int!=0){
//         printf("wall_interaction           %lf   if >0 attraction, <0 repulsion; scales equal to gravity\n", pot.wall_int);
//     }
   
//     printf("\nParticle Properties:\n");
//     print_particle_properties();


   
//     printf("\nAnalysis\n");
//     printf("bond_breakage       %.2lf\n", analysis.bond_breakage);
//     printf("RDF                 %d\n", analysis.rdfanalysis);
//     printf("xy_print            %d\n", analysis.xy_print);
//     printf("print_trajectory    %d\n", analysis.print_trajectory);
//     printf("bond_tracking       %d\n", analysis.bond_tracking);

    
//     printf("s_histogram         %d\n", analysis.s_distribution);
//     printf("cluster_analysis    %d\n", cluster.analysis);
//     printf("adjacency printing  %d\n",analysis.adjacency); 
//     printf("bond average        %lf\n", analysis.bond_avg );
//     printf("bond def. E         %lf\n", pot.bond_cutoffE);

//     printf("\n");

//     return;
// }

// void print_particle_properties(){
//     int s,n;
//     Particletype *part_type;

//      for(n=0;n<sys.nparticle_types;n++){
//         // each particletype also has a site type (may all be the same though)
//         // a particle (yet) cannot have different sites; i.e. all sites on 1 particles are the same

//         part_type=&sys.particletype[n];
//         printf("\n* particle_type %d \n",n);
//         printf("%d particle(s) have %d patche(s) with pw=%.2lf degrees: \n", part_type->nparticles, part_type->nsites, site[n].delta_degree);
//         printf("and patch vectors :\n" );

//         for(s=0; s<part_type->nsites;s++){
//             printf("   ( %10.5lf  %10.5lf  %10.5lf ) \n", part_type->site[s].x, part_type->site[s].y, part_type->site[s].z);
//             // vprint(part_type->site[s].r);
//         }
//         if(part_type->activity==1){
//             printf("particle has active force with size %lf:\n",part_type->F_A);
//             printf("   ( %10.5lf  %10.5lf  %10.5lf ) \n", part_type->e_A.x, part_type->e_A.y, part_type->e_A.z);
//         }


//         // here the patch switchi is "old " (from a nature paper, I will look up the citation), or its "new" which means explicit integration 
            
//         if (site[n].s_accent==0){
//             printf("S old:  ");
//             printf("Saccent = 0.5*(1.0-cos(PI*[costheta-cosdelta_particle]/[1-cosdelta_particle])) \n");
//             printf("        = 0.5*(1.0-cos(PI*[costheta-%10.5lf]/%10.5lf)) \n",site[n].cosdelta, site[n].oneover_cosdelta);
//         }
//         else if (site[n].s_accent==1){
//             //see  APPENDIX B.3 of J. Chem. Phys. 155, 034902 (2021); doi: 10.1063/5.0055012  
//             printf("S integrated:   ");
//             printf("Saccent = exp[-cx^2 -dx^3 -ex^4 -fx^5 -gx^6 -hx^7 -ix^8] \n");
//             printf("  c                 %.18lf\n", site[n].s_int.switch_expc);
//             printf("  d                 %.18lf\n", site[n].s_int.switch_expd);
//             printf("  e                 %.18lf\n", site[n].s_int.switch_expe);
//             printf("  f                 %.18lf\n", site[n].s_int.switch_expf);
//             printf("  g                 %.18lf\n", site[n].s_int.switch_expg);
//             printf("  h                 %.18lf\n", site[n].s_int.switch_exph);
//             printf("  i                 %.18lf\n", site[n].s_int.switch_expi);
//             if (site[n].s_int.switch_unit==0){
//                 printf("x in degrees\n");
//             }
//             else if (site[n].s_int.switch_unit==1){
//                 printf("x in cosangle\n");
//             }
//         }
//         else if (site[n].s_accent==2){
//             printf("S linear:   ");
//             printf("Saccent = 1-x/theta_p \n");
//             printf("Saccent = 1-x/%5.3lf \n",site[n].delta_degree);
//         }
//         else if (site[n].s_accent==3 && (site[n].S_fixed>0 && site[n].S_fixed<=1)){
//             // you can only make the combination of S-accent=3 and S_fixed > 0 and <=1
//             printf("S' is fixed and is set to %.2lf\n", site[n].S_fixed );
//         }
//         printf("the threshold of calculating S:\n");
//         printf("  cosdelta          %lf  \n", site[n].cosdelta);
//         printf("  delta             %lf  \n", site[n].delta_degree);
//         if(sys.gravity>0){
//             printf("and gravity parameters:\n");
//             printf("  gravity factor    %lf \n", sys.gravity);
//             printf("  zcut              %lf [sigma]\n", part_type->zcut);
//             printf("  radius            %lf [sigma]\n", part_type->radius);
//             printf("  Fg                %lf [[kT/sigma]]\n", part_type->fg);
//             printf("  b_zc              %lf [sigma]\n", part_type->b_zc);
//         }
        
//     }
//     printf("switch_method      %d (0 is S0, 1 is S90, 2 is lincomb, 3 is linear)\n\n", sys.switch_method);
    
//    return;
// }




