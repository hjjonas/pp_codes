#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "path.h"

/*SIMONS POTENTIAL PARAMETERS*/
//Vrep
double potential_repulsive_energy_sdist(double s){
    /* this is for R2*/
    double erep;


    erep= pot.Arep *exp(-s/pot.dl);
    return erep;
}
//Frep
double bond_repulsive_force(double s){
    // the repulsive force of the bond 
    double Vrep, Frep;

    if (s< pot.s_forcecut){
        // F(s) = F(s_forcecut) + F'(s_forcecut)*(s_forcecut-s)
        Frep = potential_repulsive_energy_sdist(pot.s_forcecut)/pot.dl + second_der_Vrep(pot.s_forcecut)*(pot.s_forcecut-s);
    }
    else{
        Frep=potential_repulsive_energy_sdist( s)/pot.dl;

    }

    return Frep;
}
//d^2Vrep/ds^2
double second_der_Vrep(double s){
    /*2nd derivative of Vrep ->> d^2Vrep/ds^2*/
    double Vrep,snd_der;

    Vrep=potential_repulsive_energy_sdist( s);
    snd_der = Vrep/(pot.dl*pot.dl);

    return snd_der;
}


//Vc
double potential_attractive_energy_sdist(double s){
    /* this is for R2*/
    double Uatr;

    Uatr = - (pot.A/pot.xi) * exp(-(s*s)/(pot.xi*pot.xi));

    return Uatr;
}

//Fc
double bond_attractive_force(double s){
    // the attractive force of the bond 
    double Vc;
    Vc=potential_attractive_energy_sdist( s);

    return 2.*s/(pot.xi*pot.xi)*Vc;
}

//d^2Vc/ds^2
double second_der_Vc(double s){
    /*2nd derivative of Vc ->> d^2Vc/ds^2*/
    double Fc, Vc,snd_der, s2_xi2= (s*s)/(pot.xi*pot.xi);

    Fc=bond_attractive_force(s);
    Vc=potential_attractive_energy_sdist( s);
    // snd_der = -2.*Vc/(pot.xi*pot.xi*pot.xi) + (2.*s)/pot.xi*(1.*Fc);
    snd_der = Vc * (4.*s2_xi2-2./(pot.xi*pot.xi));

    return snd_der;
}

void setup_Simons_potential_parameters(void){
    /*determine here A and xi. You choose :
    choose dt only as:  0.12 upto 0.22
    choose r_wetting only as:  0.40 upto 0.56
    choose surface_charge only as:  -0.10 upto -0.38

    these parameters are based on a system with:
    salt concentration:             1mM MgSO4   = 4.0 
    volume fractie lutidine =       25%vol      
    DP-A particles met diameter=    3.2 micron 
    patch curvature R=              1.0 micron      
        */
    double dt,rwet,sc;
    double A,xi,Arep,Apart1,Apart2,xipart1,xipart2;

    pot.dl=0.0008673892526; // = 2.775645608174734e-03 micron/ 3.2 micron

    printf("    SIMONS POTENTIAL PARAMETERS:\nthese parameters are based on a system with:\n");
    printf("salt concentration:             1mM MgSO4   = 4.0 \n");
    printf("volume fractie lutidine =       25\%%vol      \n");
    printf("DP-A particles met diameter=    3.2 micron \n");
    printf("patch curvature R=              1.0 micron  \n");
    printf("Debye length    =               2.775645608174734e-03 micron / %20.15lf sigma \n",pot.dl);


    printf("these parameters are fitted between:\n");
    printf("rwet in [0.30,0.60]\n");
    printf("sc   in [-0.50,-0.01] \n");
    printf("dt   in [0.10,0.22]\n");

    // so read the input parameters, 
    dt=pot.dT;
    rwet=pot.r_wetting;
    sc=pot.surface_charge;

    if(dt<0.10 || dt>0.22){
        printf("pot.dT= %lf  \n", dt);
        error(" choose only dt between 0.10 and 0.22  ");
    }
    if(rwet<0.30 || rwet>0.60){
        printf("pot.r_wetting= %lf  \n", rwet);
        error(" choose only rwet between 0.30 and 0.60");
    }
    if(sc>-0.01 || sc<-0.50){
        gprint(pot.surface_charge);
        error("choose surface_charge between -0.50 and -0.01\n");
    }
    
    


    /*from files: 
    fit_A_poly_4x4dt_patchrad1.00mum
    fit_xi_poly_4x4dt_patchrad1.00mum.txt
    fit_Arep_poly_3x0dt_patchrad1.00mum.txt*/

    /* pot.sigma_mum=1.0; //[mum]
    Apart1=0.056782502496567314  +rwet*0.260215676538175344  +rwet*rwet*-2.087133590729501886  +rwet*rwet*rwet*14.759587560733907097  ;    
    Apart2=0.011502266054040660  +dt*5.278593996230584118  +dt*dt*-29.955521677043037698  +dt*dt*dt*51.198280765408576087  ;
    pot.A=Apart1*Apart2;//[dimensionless]

    xipart1=0.004056363352597567  +rwet*-0.079426720163945003  +rwet*rwet*0.031859894515664684  +rwet*rwet*rwet*-0.011776312522848912  ;   
    xipart2=-0.486451801907159032  +dt*1.964137077233496065  +dt*dt*-3.757782035401259879  +dt*dt*dt*0.711773304860235934 ;
    pot.xi=xipart2*xipart1; //   [dimensionless]
    
    pot.Arep=0.000084830341705028  +sc*0.002376164488813220  +sc*sc*1295462.677729546325281262    ;*/

    /*in directory update_II_200826_1.0mum*/
    Apart1=-0.015061862283229362  +rwet*2.179971050601182458  +rwet*rwet*-9.919260032900107049  +rwet*rwet*rwet*51.479504724036686980  ;
    Apart2=-0.226173154804393844  +dt*7.246060901785060793  +dt*dt*-60.600308320350059432  +dt*dt*dt*221.602339052154405863  +dt*dt*dt*dt*-304.832948001714555630  ;

    pot.A=Apart1*Apart2;// [dimensionless]

    xipart1=0.219318405905894714  +rwet*-6.415286942580221030  +rwet*rwet*0.684058978296809062  +rwet*rwet*rwet*0.685914887127846429  ;
    xipart2=-0.007218874637425805  +dt*0.071211357421138557  +dt*dt*-0.501629147203333958  +dt*dt*dt*1.870137054939489119  +dt*dt*dt*dt*-2.783413465699470013  ;

    pot.xi=xipart2*xipart1;

    pot.Arep=-0.000060954325068766  +sc*-0.000689805950692273  +sc*sc*1295484.309107528068125248  ;



    gprint(pot.A);
    gprint(pot.xi);
    gprint(pot.Arep);
    return;


}
