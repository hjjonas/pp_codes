#include "path.h"


double S_value(double cositheta, double cosjtheta,double cos_ijtheta, int ptypei, int ptypej){
    /* the definition of S and S' define the effective switch funciton S;
    S =max(S'i,S'j) = S0        option 0 of site[ptype].switch_method
    S =S'i*S'j  = S90           option 1
    S = lincomb SO, S90, S180   option 2 
    
    S' is smooth switch         option 0 of site[ptype].s_accent;
    S' is integrated switch     option 1
    S' is linear function       option 2*/
    // returns the value of S given the 3 angles between two particles 

    // see J. Chem. Phys. 155, 034902 (2021); doi: 10.1063/5.0055012 for details
    double Si, Sj, S;
    

    if ((ptypei>sys.nparticle_types) || (ptypej>sys.nparticle_types)){
        dprint(ptypei);
        dprint(ptypej);
        dprint(sys.nparticle_types);
        error("S_value the Ptype is too large");
    }

    double cosdelta_i=site[ptypei].cosdelta, cosdelta_j=site[ptypej].cosdelta;

    if ( (cosdelta_i>cositheta) || (cosdelta_j>cosjtheta) ){
        return S=0;
    }
    
    // then create S
    switch (sys.switch_method) {
        case 0:
            // first calc Si' and S'j
            Si = Saccent(cositheta, ptypei);
            Sj = Saccent(cosjtheta, ptypej);
            S = S0(Si, Sj);
            break;
        case 1:
            // first calc Si' and S'j
            Si = Saccent(cositheta, ptypei);
            Sj = Saccent(cosjtheta, ptypej);
            S = S90(Si, Sj);
            break;
        case 2:
            /*the linear combination method*/
            S = switch_method_2(cositheta, cosjtheta, cos_ijtheta, ptypei, ptypej);
            break;
        default:
            error("choose switch_method = 0 (S0), 1(S90), 2(lincomb)");
    }

    if((S<-1e-1) || (S-1.>1e-3 && S>1.) || ((S/S!=1) && S>0) ){
        printf("WARNING: S is either S<-0.10 or S>1.001 or S=Nan\n");
        printf(" Si        = %20.18lf\n", Si);
        printf(" Sj        = %20.18lf\n", Sj);
        printf(" S         = %20.18lf\n", S);
        printf(" cositheta = %20.18lf\n", cositheta);
        printf(" cosjtheta = %20.18lf\n", cosjtheta);

        printf(" S0(Si,Sj) = %20.18lf\n", S0(Si,Sj));
        printf(" S90(Si,Sj) = %20.18lf\n", S90(Si,Sj));
        printf(" S180(Si,Sj) = %20.18lf\n", S180(Si,Sj));

    }
    return S;
}

double Saccent(double costheta, int ptype ){
    //returns the S'(i) of paritcle i with angle costheta and particletype ptype
    double Saccent;

    switch (site[ptype].s_accent) {

        case 0:
            // printf("CALC smooth_S\n");
            Saccent = smooth_S(costheta, ptype);
            break;
        case 1:
            // printf("CALC integrated_S\n");
            Saccent = integrated_S(costheta, ptype);
            break;
        case 2:
            // printf("CALC linear S\n");
            Saccent = linear_S(costheta, ptype);
            break;
        case 3:
            Saccent = site[ptype].S_fixed;
            break;
        default:
            gprint(costheta);
            dprint(ptype);
            dprint(site[ptype].s_accent);
            error("particletype[ptype].s_accent can only be 0 (for Ssmooth) and 1(for Sintegrate) and 2 for linear S, and 3 for fixed");

        
    }   
    return Saccent;

}

double linear_S(double costheta, int ptype){
    double theta_p=site[ptype].delta_degree,  Slin, angle;
    
    angle=(costheta<1.)? acos(costheta)*180./PI:0;
    Slin = 1.-angle/theta_p;
    return Slin;
}


double smooth_S(double costheta, int ptype){
    double Ssmooth, cosdelta_particle = site[ptype].cosdelta,oneover_cosdelta = site[ptype].oneover_cosdelta;

    Ssmooth=0.5*(1.0-cos(PI*(costheta-cosdelta_particle)*oneover_cosdelta));
    return Ssmooth;
}

double integrated_S(double cosangle, int ptype){
    /* switch_expbcd is a switch function as integrated with the patch integration scheme,
    and fitted to a curve exp(cx^2 + dx^3 + ex^4 + fx^5 + gx^6) where x =1-cos(theta) and theta the angle
    the first parameter in the exponent is c, because, a=0 (such that f(0)=1) and b=0, because f'(0)=0
    this is to make the switching function continuous in the fist derivative */
    double angle1, angle2, angle4;
    double cx2, dx3, ex4, fx5,gx6,hx7,ix8;
    double a, S;

    if(site[ptype].s_int.switch_unit==0){
        //use degree
        if (cosangle<=-1){
            angle1=0;
        }
        angle1=(cosangle<1.)? acos(cosangle)*180./PI:0;
    }
    else if (site[ptype].s_int.switch_unit==1){
        //use cosangel
        angle1=1.-cosangle;
        
    }
    angle2=angle1*angle1;
    angle4=angle2*angle2;

    cx2=angle2* site[ptype].s_int.switch_expc;
    dx3=angle2*angle1* site[ptype].s_int.switch_expd;
    ex4=angle4* site[ptype].s_int.switch_expe;
    fx5=angle4*angle1* site[ptype].s_int.switch_expf;
    gx6=angle4*angle2* site[ptype].s_int.switch_expg;
    hx7=angle4*angle2*angle1* site[ptype].s_int.switch_exph;
    ix8=angle4*angle4* site[ptype].s_int.switch_expi;
    a=0.-cx2-dx3-ex4-fx5-gx6-hx7-ix8;
    S=exp(a);
    if (S>0 & S/S!=1){
        gprint(a);
        gprint(cx2);
        gprint(dx3);
        gprint(ex4);
        gprint(fx5);
        gprint(gx6);
        gprint(hx7);
        gprint(ix8);
        gprint(cosangle);
        gprint(angle1);
        dprint(site[ptype].s_int.switch_unit);
    }

    return S;
}

double S0(double Si,double Sj){
    //kleinste waarde voor S bepaalt
    return min(Si,Sj);
}

double S90(double Si,double Sj){
    return Si*Sj;
}

double S180(double Si,double Sj){
    double S;
   
    S=S90(Si,Sj) +0.7*(S90(Si,Sj)-S0(Si,Sj));
    return S;
}

double switch_method_2(double cositheta, double cosjtheta, double cos_ijtheta, int typei, int typej){
    /*the linear combination method*/
    double angle1, angle2, angle3,S,costheta3;
    double lambda1, lambda2;
    double Si,Sj;
    // include theta3 for a linear combination of cube S0, sphere S90, and the combi S180.
        /*theta3 [0,180], subtract projection of rij of p1 and p2 and take vector_inp*/
    
    
    Si = Saccent( cositheta, typei);
    Sj = Saccent( cosjtheta,  typej);

    if(cositheta==1.0 || cosjtheta==1.0){
        // if one of the angles is zero, the rotation around angle3 is invariant. Thus, angle3 is arbitraty 
        S=max(Si,Sj);
    }
    else{
        // error("do not use switch method3")

        /* angle3=nan if costheta3>1 or costheta3<-1 (can sometimes happen)*/
        if((cos_ijtheta>=1) || (cos_ijtheta<=-1)){
                if(cos_ijtheta>=1.) angle3=0.;
                if(cos_ijtheta<=-1.) angle3=180.;
        }
        else{   angle3=acos(cos_ijtheta)*180./PI;
        }

        
        if(angle3>=0. && angle3<=90.0){
            lambda1=angle3/90.;
            S=lambda1*S90(Si,Sj)+(1.-lambda1)*S0(Si,Sj);

        }
        else if(angle3>90.0 && angle3<=180.0){
        // else if(cos_ijtheta<0 && cos_ijtheta>=-1){
            lambda2=(angle3-90.)/90.;
            S=(1.-lambda2)*S90(Si,Sj)+lambda2*S180(Si,Sj);
        }
        else{
            printf("****esle: \n");
            printf("cositheta = %.18lf \n", cositheta);
            printf("cosjtheta = %.18lf \n", cosjtheta);
            printf("cos_ijtheta = %.18lf \n", cos_ijtheta);
            gprint(angle3);

            error("angle3 is bigger than 180??");
        }
    }
    return S;
}

double dS_dcostheta(Slice *psl, double cositheta, int ipart, double cosjtheta, int jpart ){
    /*write the dS/dcostheta there for various switch function 
    take the derivative to cositheta

    S(i,j) or S(i) something else
    dS/dcosi = S(j)* dS(i)/dcosi*/

    double dS_dcostheta;
    int ptypei, ptypej;
    ptypei=psl->pts[ipart].ptype;
    ptypej=psl->pts[jpart].ptype;


    
    // switch_method is now in site[ptype].switch_method
    if (sys.switch_method==0){
        dS_dcostheta=dS0_dcostheta(cositheta,ptypei,cosjtheta,ptypej);  
    }
    else if(sys.switch_method==1){
        dS_dcostheta=dS90_dcostheta(cositheta,ptypei,cosjtheta,ptypej);
    }
    else if(sys.switch_method==2){
        error("im not sure how the derivative of this one goes. Im will probably not use it, make it if  you want to.");
    }
    else{
        error("choose switch_method = 0,1,2");
    }


    return dS_dcostheta;
}

double dSaccent_dcosangle(double cosangle, int ptype){
    /* the derivative wrt to cosangle of Saccent 
    S'=exp(a(x))
    dS'/dcostheta   = dS'/da * da/dx * dx/dcostheta = 
                    = S' * b * dx/dcostheta 
    with a=0-cx2-dx3-ex4-fx5-gx6-hx7-ix8 and x dependent on choice of switch_varialble*/

    double angle1,angle2,angle4;
    double b,dSaccent_dcosangle,dx_dcostheta;
    double cx1_2, dx2_3, ex3_4, fx4_5, gx5_6,hx6_7,ix7_8;

    // make distinction between particletype
    if (site[ptype].s_accent==0){ // from R. Guo, J. Mao, X.-M. Xie, and L.-T. Yan, Sci. Rep. 4, 7021 (2015).
        dSaccent_dcosangle=dSsmooth_dcostheta(cosangle,ptype);
    }
    else if (site[ptype].s_accent==1){
        dSaccent_dcosangle= dSintegrated_dcostheta(cosangle, ptype);
    }
    else if (site[ptype].s_accent==2){
        error("if s_accent=2 the switchfunction is linear and discontinuous");
    }
    else if (site[ptype].s_accent==3){
        dSaccent_dcosangle=0;
    }
    

    

    return dSaccent_dcosangle;
}

double dSintegrated_dcostheta(double cosangle, int ptype){
    /* the derivative wrt to cosangle of Saccent 
    S'=exp(a(x))
    dS'/dcostheta   = dS'/da * da/dx * dx/dcostheta = 
                    = S' * b * dx/dcostheta 
    with a=0-cx2-dx3-ex4-fx5-gx6-hx7-ix8 and x dependent on choice of switch_varialble*/

    double angle1,angle2,angle4;
    double b,dSaccent_dcosangle,dx_dcostheta;
    double cx1_2, dx2_3, ex3_4, fx4_5, gx5_6,hx6_7,ix7_8;

    if(site[ptype].s_int.switch_unit==0){
        if (cosangle<=-1){
            angle1=180.;
        }
        if (cosangle>=1.){           
            angle1=0.;
            dx_dcostheta=0;
        }
        else{                        
            angle1=(cosangle<1.)? acos(cosangle)*180./PI:0;
            dx_dcostheta = -1./sin(angle1*PI/180.)*(180./PI);
        }

        /*dx/dcostheta= dx/dtheta * dtheta/dcostheta = dx/dtheta * (dcostheta/dtheta)^(-1)
        x=theta */
    }
    else if (site[ptype].s_int.switch_unit==1){
        angle1=1.-cosangle;
        //x=1-costheta so dx/dcostheta = -1
        dx_dcostheta=-1.;
        
    }

    angle2=angle1*angle1;
    angle4=angle2*angle2;

    /* b=da/dx  */
    cx1_2=2.* site[ptype].s_int.switch_expc*angle1;
    dx2_3=3.* site[ptype].s_int.switch_expd*angle2;
    ex3_4=4.* site[ptype].s_int.switch_expe*angle2*angle1;
    fx4_5=5.* site[ptype].s_int.switch_expf*angle4;
    gx5_6=6.* site[ptype].s_int.switch_expg*angle4*angle1;
    hx6_7=7.* site[ptype].s_int.switch_exph*angle4*angle2;
    ix7_8=8.* site[ptype].s_int.switch_expi*angle4*angle2*angle1;

    b = 0.-cx1_2-dx2_3-ex3_4-fx4_5-gx5_6-hx6_7-ix7_8;
    
    dSaccent_dcosangle = Saccent(cosangle, ptype) * b * dx_dcostheta;

    // gprint(Saccent(cosangle, ptype));
    // gprint(b);
    // gprint(dx_dcostheta);
    return dSaccent_dcosangle;
}

double dSsmooth_dcostheta(double costheta , int ptype){
    double dSsmooth_dcostheta;
    double cosdelta_particle = site[ptype].cosdelta,oneover_cosdelta = site[ptype].oneover_cosdelta;

    dSsmooth_dcostheta=HALFPI*sin(PI*(costheta-cosdelta_particle)*oneover_cosdelta)*oneover_cosdelta;
    // gprint(dSsmooth_dcostheta);
    // gprint(cosdelta_particle);
    // gprint(oneover_cosdelta);
    return dSsmooth_dcostheta;
}

double dS0_dcostheta(double cosangle1, int ptype1, double cosangle2, int ptype2){
    /*derivative of S0 wrt cosangle1 */
    double dS0_dcostheta;

    /*check which angle is biggest (i.e. cosangle smallest) and use that one as theta*/
    if(cosangle2<cosangle1){
        dS0_dcostheta=0;
    }
    else{
        dS0_dcostheta=dSaccent_dcosangle(cosangle1, ptype1)*Saccent(cosangle2,ptype2);
    }
    
    return dS0_dcostheta;
}

double dS180_dcostheta(double cosangle1, int ptype1, double cosangle2, int ptype2){
    double dS180_dcostheta,dS90=dS90_dcostheta(cosangle1,ptype1,cosangle2,ptype2);
    double dS0 = dS0_dcostheta(cosangle1,ptype1,cosangle2,ptype2);
   
    dS180_dcostheta=dS90 +0.7*(dS90-dS0);
    return dS180_dcostheta;
}

double dS90_dcostheta(double cosangle1, int ptype1, double cosangle2, int ptype2){
    /*derivative of S90 wrt cosangle1
    with S90 = S'(cosangle1) * S'(cosangle2)  */
    double dS90_dcostheta;

    /*dS90_dcostheta = dS'/dcostheta1 *S'(costheta2)*/
    dS90_dcostheta = dSaccent_dcosangle(cosangle1,ptype1)*Saccent(cosangle2,ptype2);
    
    return dS90_dcostheta;
}
