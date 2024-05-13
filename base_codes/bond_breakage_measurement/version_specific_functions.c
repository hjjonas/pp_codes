#include "path.h"

//################################################################################################################
//######################################### BOND BREAKAGE FUNCTIONS ##############################################
//################################################################################################################

/*------------------LOCAL FUNCTIONS------------------------------------------*/

void set_transitions_to_zero(int );
double return_Svalue(Slice *, int  , int  );
void update_histogram_S(double *, double *);
int  return_currect_loc(Slice *, int, int);
void breakage_probabilities(int );
void append_to_file(char [100], double , int );
void print_to_file(char [100], double ,int );
void printing_Svalues(int);
void print_hist(char [100], double *, int, int  );
void print_transitions_all_rings(Slice *psl);
void print_transitions_all(Slice *psl);

void set_reduced_bondnumber_chain(Slice *);
void set_reduced_bondnumber_ring(Slice *);

int check_interface_crossing(Slice * , int , int );
void bond_broken(Slice *);
/*---------------------------------------------------------------------------*/

TST breakage_all_tst[NPART],breakage_i_tst[NPART];
int bondnumbers[NPART][NPART], reduced_bondnumbers[NPART][NPART], Nred_bonds;
TST_Lambda tst_lambda;
Statistics tau_sims[NPART];
/*---------------------------------------------------------------------------*/

// use these local dummies if you want to print stuff to screen for testing/debuggin
int test_printing=  0;
int test_printing_2=0;


int bondnumbers[NPART][NPART]={-2}, reduced_bondnumbers[NPART][NPART]={-2}, Nred_bonds;
TST_Lambda tst_lambda;


void version_specific_analysis(Slice *psl){
	// this funciton is located inside terminate_block() at the outerloop in main.c	
    return ;
}

void every_timestep_analysis(Slice *psl){
    // perform these analysis function after every timestep. 
    // practically only used for bond breakage analysis or making histogrmas

    // see DOI: 10.1039/d3sm01255g
    int brokencheck, do_printing=0;
    if (langevin.bond_breakage_analysis!=1){
        return;
    }

    //loop over unique particle pairs, check if they are bond in the start configuration
    for ( int ipart=0; ipart<psl->nparts;ipart++){ 
        for ( int jpart=ipart+1; jpart<psl->nparts;jpart++){ 
            // if they made a bond in the start configuration, check the current status of that particle pair
            if (bond_check(start_slice,ipart,jpart)){
                brokencheck=check_interface_crossing(psl, ipart, jpart); // check for crossings every timestep
                if (brokencheck){ 
                    do_printing=1;
                    printf(" >> broken at bond %d = reduced bond %d \n",bondnumbers[ipart][jpart],reduced_bondnumbers[ipart][jpart] );
                }
            }
        }
    }
    // if any bond has broken, stop the simulation 
    if (do_printing) bond_broken(psl); do_printing=0;
    
}

void innerloop_analysis(Slice *psl){
	// this funciton is located inside in innerloop in main.c
	// so after ncycle2 bmd, mc or trajectory reads, you perform this function. 
    return;
}


//######################################################################################################################################
//######################################### TST: CALCULATIONS OF INTRINSIC RATES #######################################################
//######################################################################################################################################

/* 	works for a chain of npart-1 bonds, the bonds numbers should be ordered ! (use start_type =2 for making a chain)
    or when analysing the ring, you first had to identify the reduced bond number (done in init.c) 

  	so, BMD is just doing its thing by moving the particles around in the potential well

  	you need to track when you cross an interface, these interfaces are defined based on the interparticle distance (or energy? )
  	interface definitions: 		1st interface is at 		E_lambda=0.7*Emin_lambda		(lambda_12)
								2nd interface is at 		E_lambda=0.06*Emin_lambda		(lambda_23)
								3rd interface is at         E_lambda=0 & r=0.5 sigma    	(lambda_34)
    
    see for more info/definitions the Ref https://pubs.rsc.org/en/content/articlelanding/2024/sm/d3sm01255g
    DOI https://doi.org/10.1039/D3SM01255G 

	so I need a function that checks if I crossed an interface (done in check_interface_crossing), in positive or negative manner (left to right or reverse forward backward )
	if the interface didn't change always update the tst[ipart].times->last to keep track of the last time you were at a certain region (region 1 2 3 4), 
    and if the crossing in postitive direction (increasing the region number).  add it to tst[ipart].timess.count and

	if the interface changes, you need to register the interface change in tst[ipart].transition_i.N_13 ++ (or N1, N_24, N_31 etc)
	and change the flag from were you came from (which region) no next time you perform    check_interface_crossing, you know your history                     

    8 mei 2023 - What you should  measure is <k> = sum N_12 / sum time * sum N_1->3 / sum N_12 *  sum N_2->4 / sum N_1->3
    So that is what I print to a file in the function

 */

int check_interface_crossing(Slice *psl , int ipart, int jpart){
    // this function (check_interface_crossing) is performed after every timestep/bmd move!! 
    
    // given an configuration, identify in which region the configuration lies 
    int newloc=return_currect_loc(psl, ipart,  jpart);

    // for both rings and chains we use (reduced) bond numbers.
    int bond_nr=bondnumbers[ipart][jpart];

    // if both flags are zero then you just started the simulation below 0 (tst[ipart].loc=0) or between 0 and E (tst[ipart].loc=1) ) 
    // you need to start the first tst[ipart].time struct and define the flags
    if ((breakage_i_tst[bond_nr].flag_from1==0) && (breakage_i_tst[bond_nr].flag_from3==0)){
        // printf("in beginning. first go into tst[bond_nr].loc=%d . current energy=%lf (newloc=%d) <tst[bond_nr].lambda_12=%lf\n",tst[bond_nr].loc,psl->energy,newloc,tst[bond_nr].lambda_12);
        if ((breakage_i_tst[bond_nr].loc==1) && (newloc==2)){        // >>>>>>>>> positive crossing of  lambda_12  <<<<<<<<<
            printf("    START observing for bond %d  ... ",bond_nr);
            breakage_i_tst[bond_nr].flag_from1=1;       // you enter for the first time the 0 boundary (only possible at the beginning of the simulations )
            breakage_i_tst[bond_nr].last_lambda_12_time=psl->c_time; 
        }
        breakage_i_tst[bond_nr].loc=newloc;
        return 0 ;
    }

    if (newloc!=breakage_i_tst[bond_nr].loc){   // you crossed an interface, possible crossings : types
                            // 1->2 : (positive) recrossing
                            // 2->3 : (positive) recrossing or a first crossing (i.e. flag_from1 change), 
                            // 2->1 : (negative) recrossing, 
                            // 3->2 : (negative) recrossing,
                            // 3->4 :  official bond breakage
        if      ((breakage_i_tst[bond_nr].loc==1) && (newloc==2)){      // >>>>>>>>> positive crossing of  lambda_12 <<<<<<<<<
            if (breakage_i_tst[bond_nr].flag_from1){        // recrossing
                if (test_printing_2)  printf(" positive crossing of  lambda_12 ,");
                
                //  FLUX CALCULATION
                // before you update the tst[bond_nr].times.last_lambda_12_time, calculate the invflux_12
                double timedif=psl->c_time - breakage_i_tst[bond_nr].last_lambda_12_time;
                running_statistics(&breakage_i_tst[bond_nr].invflux_12,timedif); 
                     
                // update the values for next flux measurement
                breakage_i_tst[bond_nr].N_1++;  
                breakage_i_tst[bond_nr].last_lambda_12_time=psl->c_time;  
            }
            else if (breakage_i_tst[bond_nr].flag_from3){
                dprint(breakage_i_tst[bond_nr].flag_from1);
                dprint(breakage_i_tst[bond_nr].flag_from3);
                dprint(breakage_i_tst[bond_nr].loc);
                dprint(newloc);

                error("its not possible to have flag_from3 and loc change from 1->2 which means positive lambda_12 crossing");
            }          
        }
        else if ((breakage_i_tst[bond_nr].loc==2) && (newloc==3)){      // >>>>>>>>> positive crossing of  lambda_23 <<<<<<<<<
            if (breakage_i_tst[bond_nr].flag_from3){        // recrossing
                if (test_printing_2) printf(" positive crossing of lambda_23,\n");
            }
            else if (breakage_i_tst[bond_nr].flag_from1){   // region change from 1 to 3 
                if (test_printing_2) printf(" region change from 1 to 3 ,\n");

                // N_1->3
                breakage_i_tst[bond_nr].N_13++; // first crossing of lambda_23 after crossing lambda_12

                // update the histogram of S the accounts for all positive crossings of lambda_23 interface after crossing if coming from 1
                // n_1->3(S)
                double S_mech=return_Svalue(psl,ipart,jpart) ;
                breakage_i_tst[bond_nr].S_mech=S_mech; // save S_mech for later, maybe breaking, maybe rebinding
                update_histogram_S( breakage_all_tst[bond_nr].n_1to3 , &S_mech );

                // fill in new data
                breakage_i_tst[bond_nr].flag_from1=0;       // change flag
                breakage_i_tst[bond_nr].flag_from3=1;
            }
            // update last time you crossed lambda_23 forward
            breakage_i_tst[bond_nr].tau_sim=psl->c_time;         
      
            // update the histogram of S the accounts for all lambda_23 crossings
            // see Fig 8b
            double current_S=return_Svalue(psl,ipart,  jpart) ;
            update_histogram_S( breakage_all_tst[bond_nr].S_hist_lambda_23_all , &current_S); // positive crossing of lambda_23
        } 
        else if ((breakage_i_tst[bond_nr].loc==2) && (newloc==1)){      // >>>>>>>>> negative crossing of  lambda_12 <<<<<<<<<
            if (breakage_i_tst[bond_nr].flag_from1){        // you already crossed lambda_12 once before
                if (test_printing_2) printf(" negative crossing of  lambda_12,\n");
            }
            else if (breakage_i_tst[bond_nr].flag_from3){   // region  change from  2 to 1 
                // REBINDING
                if (test_printing_2) printf(" region  change from  2 to 1,\n");

                breakage_i_tst[bond_nr].N_31++;                     // add the N_31
                update_histogram_S(breakage_all_tst[bond_nr].n_3to1 , &breakage_i_tst[bond_nr].S_mech ); 

                breakage_i_tst[bond_nr].flag_from1=1;       // change flag
                breakage_i_tst[bond_nr].flag_from3=0;
            }
        }
        else if ((breakage_i_tst[bond_nr].loc==3) && (newloc==2)){      // >>>>>>>>> negative crossing of  lambda_23 <<<<<<<<<
            if (breakage_i_tst[bond_nr].flag_from3){        // you already crossed lambda_23 once before
                if (test_printing_2) printf(" negative crossing of  lambda_23,\n");
            }
            else if (breakage_i_tst[bond_nr].flag_from1){   // region  change from  2 to 3 
                if (test_printing_2) printf(" region  change from  2 to 3,  with flag form 1 \n");

                dprint(breakage_i_tst[bond_nr].flag_from1);
                dprint(breakage_i_tst[bond_nr].flag_from3);
                dprint(breakage_i_tst[bond_nr].loc);
                dprint(newloc);

                error("its not possible to have flag_from1 and loc change from 3->2 which means lambda_23 crossing");
            }
            // update the histogram of S the accounts for all lambda_23 crossings
            // see Fig 8b
            double current_S=return_Svalue(psl,ipart,  jpart) ;
            update_histogram_S( breakage_all_tst[bond_nr].S_hist_lambda_23_all , &current_S); // positive crossing of lambda_23

        }
        else if ((breakage_i_tst[bond_nr].loc==3) && (newloc==4)){      // >>>>>>>>> wehooo succesfull breakage      <<<<<<<<<
            // plus 1 for N_2to4 in breakage_i_tst
            breakage_i_tst[bond_nr].N_24++; 

            // place this one here because you only want to track the tau_sim of the broken bond. 
            // ithe <tau_sim> can be calulated with the info in the file TST_all_P_tau, becayse there breakage_i_tst[bond_nr].tau_sim is printed
            //  or via  <tau_sim>(bond_rn)= breakage_all_tst[bond_nr].tau_sim/breakage_all_tst[bond_nr].N_24           
            breakage_all_tst[bond_nr].tau_sim+=breakage_i_tst[bond_nr].tau_sim;
            
            // for tracking the breakage mechanism
            update_histogram_S(breakage_all_tst[bond_nr].n_2to4 , &breakage_i_tst[bond_nr].S_mech );

            // print P_AE, Pescape, Prebind to a file and update the averages
            breakage_probabilities(bond_nr);
            
            return 1;
        }
        else{
            // to catch errors 
            printf(" too fast loc change?(time=%lf, energy=%lf) \n",psl->c_time,total_energy(psl));
            dprint(breakage_i_tst[bond_nr].flag_from1);
            dprint(breakage_i_tst[bond_nr].flag_from3);
            dprint(breakage_i_tst[bond_nr].loc);
            dprint(newloc);
            error("stop");
        }
        breakage_i_tst[bond_nr].loc=newloc; // update the newloc to tst[ipart].loc
        if (test_printing) printf("  >> tst[%d].loc updated to %d (time=%lf, energy=%lf) \n", bond_nr, breakage_i_tst[bond_nr].loc,psl->c_time, total_energy(psl));
    }
    
    return 0;
    
}

void bond_broken(Slice *psl){
    printf("     wehooo succesfull breakage at ctime=%lf \n",psl->c_time);
    int bond_nr,redbond_nr;
    
    // ok few things to be done: 
    // printing the last time that you crossed the lambda23 interface which is bond lifetime
    
    // ok now restart from a  bound configuration    
    restart_conf(psl);

    /* print final timestamps to file and update global the escape time, transition time, and the conditional probabilities */
    // for (int ipart=0;ipart<psl->nparts-1;ipart++){
    for ( int ipart=0; ipart<psl->nparts;ipart++){ 
        for ( int jpart=ipart+1; jpart<psl->nparts;jpart++){ 
            //loop over all particle combinations, check if they are bond in start configurations
            if (bond_check(start_slice,ipart,jpart)){
                bond_nr=bondnumbers[ipart][jpart];
                redbond_nr=reduced_bondnumbers[ipart][jpart];

                // printf(" bond_nr %d\n", bond_nr);
                // printing flux throught the lambda_12 surface
                print_to_file("invflux_12", breakage_i_tst[bond_nr].invflux_12.mean, redbond_nr);

                // save to information form the single bond to all bonds 

                breakage_i_tst[bond_nr].loc=return_currect_loc(psl, ipart,  jpart); // calculate, the newloc =1
                set_transitions_to_zero(bond_nr );      //and reset the flags + counts 
                printing_Svalues(  bond_nr );

            }
        }
    }  


    if (start_slice->nbonds==start_slice->nparts){        print_transitions_all_rings(psl);}
    else if(start_slice->nbonds==start_slice->nparts-1){  print_transitions_all(psl);}

}



void printing_Svalues( int bond_nr ){
    FILE *fp;
    char filename[100];
    // stuff from all breakages and non-breakages together 
    // keeps updating the output file as long as simulation is running


    // printing (appending) the S values to file
    char extension[MAX_EXT_LENGTH],value[MAX_VALUE_LENGTH];
    sprintf(extension, "_bond%d",bond_nr);
    
    // print mechanism histograms to file. used for e.g. Fig 9

    print_hist("S_hist_n_2to4", breakage_all_tst[bond_nr].n_2to4, NBINS_MECH, bond_nr );
    print_hist("S_hist_n_1to3", breakage_all_tst[bond_nr].n_1to3, NBINS_MECH, bond_nr );
    print_hist("S_hist_n_3to1", breakage_all_tst[bond_nr].n_3to1, NBINS_MECH, bond_nr );
    print_hist("S_hist_lambda_23", breakage_all_tst[bond_nr].S_hist_lambda_23_all,NBINS_MECH, bond_nr ); // Fig 8b
    
    // for eq 14 - 16 , resp.
    double P_1to3[NPART],P_2to4[NPART],P_mech[NPART];
    for (int s=0; s<NBINS_MECH ; s++ ){
        P_1to3[s]=breakage_all_tst[bond_nr].n_1to3[s]/breakage_all_tst[bond_nr].N_13;
        P_2to4[s]=breakage_all_tst[bond_nr].n_2to4[s]/breakage_all_tst[bond_nr].n_1to3[s];
        P_mech[s]=breakage_all_tst[bond_nr].n_2to4[s]/breakage_all_tst[bond_nr].N_24;
    }

    print_hist("P_1to3", P_1to3, NBINS_MECH, bond_nr );
    print_hist("P_2to4", P_2to4, NBINS_MECH, bond_nr );
    print_hist("P_mech", P_mech, NBINS_MECH, bond_nr );

}

void print_hist(char filename[100], double *value, int length, int ipart){
    FILE *fp;

    char filename_full[150];
    char* extension = ".out";
    
    // the full filname will be composed of 2 parts: filename+"_bondX.out"
    snprintf(   filename_full, sizeof( filename_full ), "%s_bond%d%s", filename,ipart, extension );

    // printing the S values to file
    if((fp=fopen(filename_full,"w"))==NULL) {; // 
        printf("Warning: not able to append to %s\n", filename_full);
        error("stop");
    }
    else { 
        for (int i=0;i<length ; i++){
           fprintf(fp,"%.5lf\n",value[i]); 
        }
    }
                 
    fclose(fp);
    return;
}

void append_to_file(char filename[100], double value, int ipart){
    FILE *fp;
    char filename_full[150];
    char* extension = ".out";
    
    // the full filname will be composed of 2 parts: filename+"_bondX.out"
    snprintf(   filename_full, sizeof( filename_full ), "%s_bond%d%s", filename,ipart, extension );

    if((fp=fopen(filename_full,"a"))==NULL) {; // 
        printf("Warning: not able to append to %s\n", filename_full);
        error("stop");
    }
    else { fprintf(fp,"%.10lf\n",value); }
    fclose(fp);

    return;
}

void print_to_file(char filename[100], double value, int ipart){
    // overwrite the file
    FILE *fp;

    char filename_full[150];
    char* extension = ".out";
    
    // the full filname will be composed of 2 parts: filename+"_bondX.out"
    snprintf(   filename_full, sizeof( filename_full ), "%s_bond%d%s", filename,ipart, extension );


    // printing the S values to file
    if((fp=fopen(filename_full,"w"))==NULL) {; // 
        printf("Warning: not able to write to %s\n", filename_full);
        error("stop");
    }
    else { fprintf(fp,"%.10lf\n",value); }
    fclose(fp);

    return;
}


void breakage_probabilities(int bond_nr){
	// ok if you are here, there has been a succesful breakage taken place 
	// calculate the P_lambda23_lambda12 (the conditional probability P(lambda_23 | lambda_12), Psep, Prebind, tau_AfirstE, tau_AlastE. see readme file for definitions
	double P_lambda23_lambda12, Psep, Prebind;

	// note the counts are integers
	P_lambda23_lambda12  =(double)((double) (breakage_i_tst[bond_nr].N_13)/breakage_i_tst[bond_nr].N_1) ;  // eq 12 
	Psep =(double)(breakage_i_tst[bond_nr].N_24)/breakage_i_tst[bond_nr].N_13 ;  // eq 13
	Prebind =(double)(breakage_i_tst[bond_nr].N_31)/breakage_i_tst[bond_nr].N_13 ;  // Prebind= 1-Psep

	running_statistics(&breakage_all_tst[bond_nr].P_lambda23_lambda12,P_lambda23_lambda12);
	running_statistics(&breakage_all_tst[bond_nr].Psep,Psep);
    running_statistics(&breakage_all_tst[bond_nr].Prebind,Prebind);
    running_statistics(&breakage_all_tst[bond_nr].invflux_12,breakage_i_tst[bond_nr].invflux_12.mean);

	// printing everything to file
    static int print_header=1;
    char header[MAX_VALUE_LENGTH], value[MAX_VALUE_LENGTH], filename[MAX_VALUE_LENGTH];
    
    // used for e.g. Fig. 13 and 17
    // print info about the single breakage of bond_nr
    sprintf(  filename, "TST_all_P_tau.out" );
    sprintf(value,"%.15lf,%.15lf,%.15lf,%.15lf,%.15lf,%d",P_lambda23_lambda12,Psep,Prebind,breakage_i_tst[bond_nr].invflux_12.mean, breakage_i_tst[bond_nr].tau_sim,bond_nr);
    
    if (print_header){
        sprintf(header,"#P_lambda23_lambda12,Psep,Prebind,invflux_12,tau_sim,bond");
        write_append_to_file(filename, "", 'w', header);
        print_header=0;
    }

    
	write_append_to_file(filename, "", 'a', value);

    return;
}

void print_transitions_all_rings(Slice *psl){

    // printing the sum of the number of crossings to file
    //
    FILE *fp;
    char filename[100];
    snprintf(   filename, sizeof( filename ), "transition_sums_N.out" );
    int header=1, redbond_nr,bond_nr;

    int N_1[NPART]={0},N_13[NPART]={0},N_24[NPART]={0};
    double time[NPART]={0};
    for ( int ipart=0; ipart<psl->nparts;ipart++){ 
        for ( int jpart=ipart+1; jpart<psl->nparts;jpart++){ 
            //loop over all particle combinations, check if they are bond in start configurations
            if (bond_check(start_slice,ipart,jpart)){
                
                redbond_nr=reduced_bondnumbers[ipart][jpart];
                bond_nr=bondnumbers[ipart][jpart];

                N_1[redbond_nr]+=breakage_all_tst[bond_nr].N_1;
                N_13[redbond_nr]+=breakage_all_tst[bond_nr].N_13;
                N_24[redbond_nr]+=breakage_all_tst[bond_nr].N_24;
                time[redbond_nr]+=breakage_all_tst[bond_nr].tau_sim; 
                
            }
        }
    }


    if((fp=fopen(filename,"w"))==NULL) { // overwrite/update after each simulation
        error("Warning: not able append to transition_sums_N.out\n");
    }
    else {
        for ( int redbond_nr=0; redbond_nr<Nred_bonds;redbond_nr++){ 
           
            if (header==1){
                fprintf(fp,"bond,N_1,N_1->3,N_2->4,time\n");  
                header=0;
            }

            fprintf(fp,"%d,%d,%d,%d,%.15lf\n",redbond_nr, 
                                              N_1[redbond_nr] ,
                                              N_13[redbond_nr],
                                              N_24[redbond_nr],
                                              time[redbond_nr]);   
        }
    }
    fclose(fp);

    return;
}

void print_transitions_all(Slice *psl){

    // printing the sum of the number of crossings to file
    //
    FILE *fp;
    char filename[100];
    snprintf(   filename, sizeof( filename ), "transition_sums_N.out" );
    int header=1, redbond_nr,bond_nr;

    // This is for printing the reduced bond informaion
    int N_1[NPART]={0},N_13[NPART]={0},N_24[NPART]={0};

    for ( int ipart=0; ipart<psl->nparts;ipart++){ 
        for ( int jpart=ipart+1; jpart<psl->nparts;jpart++){ 
            //loop over all particle combinations, check if they are bond in start configurations
            if (bond_check(start_slice,ipart,jpart)){
                
                redbond_nr=reduced_bondnumbers[ipart][jpart];
                bond_nr=bondnumbers[ipart][jpart];

                N_1[redbond_nr]+=breakage_all_tst[bond_nr].N_1;
                N_13[redbond_nr]+=breakage_all_tst[bond_nr].N_13;
                N_24[redbond_nr]+=breakage_all_tst[bond_nr].N_24;
                running_statistics( &tau_sims[redbond_nr], breakage_all_tst[bond_nr].tau_sim) ; 
                
            }
        }
    }   


    if((fp=fopen(filename,"w"))==NULL) { // overwrite/update after each simulation
        error("Warning: not able append to transition_sums_N_reduced_bonds.out\n");
    }
    else {
        for ( int redbond_nr=0; redbond_nr<Nred_bonds;redbond_nr++){ 
           
            if (header==1){
                fprintf(fp,"bond,N_1,N_1->3,N_2->4,<tau_sim>\n");  
                header=0;
            }

            fprintf(fp,"%d,%d,%d,%d,%lf\n",redbond_nr, 
                                              N_1[redbond_nr] ,
                                              N_13[redbond_nr],
                                              N_24[redbond_nr],
                                              tau_sims[redbond_nr].mean);   
        }
    }
    fclose(fp);

    return;
}



void update_histogram_S( double *S_hist , double *S ){
    /* this function is to update the S-histogram (with NBINS_MECH bins), that tracks the value at 
    do this if the bond has broken (you reached region 4)*/

    int bin;
    double S_bin=1./(double)NBINS_MECH;

    bin=(int)(*S/S_bin);
    if( bin==NBINS_MECH) bin=NBINS_MECH-1; // if S==1 (pretty extrodinaty, put it to bin 9) so last bin is \in [0.9,1.0]
    if ((bin>=0 )&& (bin<NBINS_MECH)){ 
        S_hist[bin]++;
        }
    else {
        gprint(*S); 
        error("you are outside the bins range of tst[ipart].S_histogram ");
    }   

    *S=-1; // I used this as a way to track possible errors 
    return;
}



void set_transitions_to_zero(int bond_nr){
    // after breakage (reaching region 4), first add them to transitions_all , then reset all of these values and start again

    // overall number of transitiosn / crossings through the boundaries
    breakage_all_tst[bond_nr].N_13 +=breakage_i_tst[bond_nr].N_13;
    breakage_all_tst[bond_nr].N_31 +=breakage_i_tst[bond_nr].N_31;
    breakage_all_tst[bond_nr].N_24 +=breakage_i_tst[bond_nr].N_24;
    breakage_all_tst[bond_nr].N_1  +=breakage_i_tst[bond_nr].N_1;
    

    // transition_i is per simulation
    breakage_i_tst[bond_nr].flag_from1=0;
    breakage_i_tst[bond_nr].flag_from3=0;

    breakage_i_tst[bond_nr].N_13=0;
    breakage_i_tst[bond_nr].N_31=0;
    breakage_i_tst[bond_nr].N_24=0;
    breakage_i_tst[bond_nr].N_1=0;

    breakage_i_tst[bond_nr].last_lambda_12_time=0;
    breakage_i_tst[bond_nr].tau_sim=0; 


    return;

}


void set_reduced_bondnumber_ring(Slice *psl){
    /* do this at beginning of simulation, it walks over the ring, 
    and gives the unique bondnumber to each bond in the ring. 
    I think a big array (staring with -1's) of npart x npart (like an adjacency matrix)
    then you fill in your particle numbers in the array, and it gives back the reduced bondnumber
    */

    //some checks first 
    if ((sys.particletype[0].nsites!=3 )&& (sys.particletype[1].nsites!=2)){
        error("first type should be TPP, and second DP particle");
    }

    // the first NTTP particles are TPP particles, followed by NDP DP particles
    // [ TPP1, TPP2, ... ,TPPn , DP1, ... , DPn]
    // the chains are DPn/TPPn long 
    // maybe first walk over the chains and fill in 
    int NDP,Lchain,Nchains,NTPP;
    NDP=sys.particletype[1].nparticles;
    NTPP=sys.particletype[0].nparticles; // first type is always TPP

    Lchain=NDP/NTPP;
    Nchains=NTPP;
    int DP1=NTPP;// the first DP particle has number NTPP
    
    // put DP particles in a list
    IntArray DP_list; 

    initIntArray(&DP_list,  (size_t) NDP);
    for(int i=NTPP; i<(NDP+NTPP); i++){        
        insertIntArray(&DP_list, i);
    }

    psl->energy=total_energy(psl);

    int DPc=DP1;
    Nred_bonds=Lchain+1;
    int bond_nr=0,redbond_nr;

    for (int chain=0; chain<Nchains; chain++ ){
        // restart; set bond number to zero   
        redbond_nr=0;
        // dprint(redbond_nr)
        // walk over chain of length Lchain
        for (int L=0; L<Lchain; L++){
            // printf(" the %d th particle = %d in chain %d of length %d\n",L,DPc,chain,Lchain );       
            // the bonds of DPc
            for (int n=0; n<psl->pts[DPc].nbonds;n++){

                int bound_particle=psl->pts[DPc].bonds[n];
                // printf("     particle %d and %d make bond  \n",DPc ,bound_particle);  
                if ((DPc<bound_particle) || (bound_particle<NTPP) ) {
                    bondnumbers[DPc][bound_particle]=bond_nr;
                    bondnumbers[bound_particle][DPc]=bond_nr;
                    bond_nr++;
                }     

                // check if it as a DP particle int the list, or a TPP particle
                if ((checkElementXIntArray(&DP_list, bound_particle )==1) || (bound_particle<NTPP) ){
                    // dprint(checkElementXIntArray(&DP_list, bound_particle ));
                    printf("      bond between (%d-%d) , given reduced bondnumebr:  %d \n",DPc,bound_particle,redbond_nr);       
                    reduced_bondnumbers[DPc][bound_particle]=redbond_nr;
                    reduced_bondnumbers[bound_particle][DPc]=redbond_nr;

                    if (bound_particle>= DP1){removeElementXIntArray( &DP_list ,  bound_particle);} 
                    redbond_nr+=1;
                }
            }
            if (checkElementXIntArray(&DP_list, DPc )==1){removeElementXIntArray( &DP_list ,  DPc);}
            DPc+=1; // go to next particle 
        }
    }

    freeIntArray(&DP_list) ;

    if (redbond_nr!=Nred_bonds) error("problem with finding reduced bond numbers for rings");

    return;
}

void set_reduced_bondnumber_chain(Slice *psl){
    int red_bondnr=0;
    
    // loop over the particle pairs, if they make a bond and 
    // note that for chains, i'm actually not using reduced bond numbers, 
    // because its easy to get a lot of data for each bond.
    
    for (int i=0; i<psl->nparts;i++){
        for (int j=i+1; j<psl->nparts;j++){
            if (bond_check(psl,  i,  j)){
                reduced_bondnumbers[i][j]=red_bondnr;
                reduced_bondnumbers[j][i]=red_bondnr;

                bondnumbers[i][j]=red_bondnr;
                bondnumbers[j][i]=red_bondnr;

                printf("  for particle %d and %d, they have bond number %d \n",i,j,red_bondnr);

                red_bondnr++;
            }
        }
    }

    Nred_bonds=red_bondnr;
    if (Nred_bonds==0 ) error("no bonds found");
    if (red_bondnr!=psl->nparts-1) error("problem with finding reduced bond numbers for chain");

    // dprint(Nred_bonds);
    return;
}


int  return_currect_loc(Slice *psl, int ipart, int jpart){
    /*  the energy/distance-landscape into regions is defined by lambda(r,S)
     */

    double r_lambda,S,r; 

    // all this for S_value
    vector rij,u1;
    double  cositheta,  cosjtheta, cos_ijtheta;
    particles_distance_length_vector(psl,ipart,jpart,&r_lambda,&r,&rij);

    int bond_nr=bondnumbers[ipart][jpart];

     // if broken, i.e., r_lambda=>sys.bondbreakage_treshold. no need to calculate energies , directly return loc
    if (r_lambda>tst_lambda.bondbreakage_treshold) { 
        return 4;
    }

    // if distance pot.s_cutoff<r_lambda<ys.bondbreakage_treshold , directly return loc; this line can also be removed, but is saves some calculations
    if (r_lambda>pot.s_cutoff){  
        return 3;
    }
    scalar_divide(rij,r,u1);
    orientation_parameters( psl,  ipart,  jpart, u1, &cositheta, &cosjtheta, &cos_ijtheta);
    S=S_value( cositheta,  cosjtheta, cos_ijtheta, psl->pts[ipart].ptype, psl->pts[jpart].ptype);

    // if distance smaller than s_min, set r_lambda to s_min 
    MAX(r_lambda,pot.s_min, r_lambda); 
    
    double E_lambda,Erep,Ec;
    // Erep=potential_repulsive_energy_sdist(r_lambda);
    Ec=potential_attractive_energy_sdist(r_lambda);
    E_lambda=Ec*S; // eq 9 of paper


    if (E_lambda<=tst_lambda.lambda_12){     
        return 1;
    }
    else if ((E_lambda>tst_lambda.lambda_12)  &&  (E_lambda<=tst_lambda.lambda_23)){   
        return 2;
    }
    else if ((E_lambda>=tst_lambda.lambda_23) ){  // in the positive region
        return 3;      
    }
    else{
        printf("\nsomething is wrong\n");
        gprint(tst_lambda.lambda_23);
        gprint(tst_lambda.lambda_12);
        gprint(E_lambda);
        gprint(Erep);
        gprint(Ec);
        gprint(S);
        gprint(r_lambda);
    }
    error("was not able to locate the current loc? ");

    return -1;
}




/*_________________________MEMORY ALLOCATION AND FREEING________________*/
void allocate_memory(Slice *psl){
	// allocate versions specific memory
	// maybe put this in init?

	// dynamically allocate memory for cluster.size_distribution.bin 
	cluster.size_histogram.bin = (Statistics *)calloc(psl->nparts, sizeof(Statistics));
	cluster.size_distribution.bin = (Statistics *)calloc(psl->nparts, sizeof(Statistics));

	cluster.size_histogram.length=psl->nparts;
	cluster.size_distribution.length=psl->nparts;


	return;
}

void free_all_memory(void){
	// free the memory you allocated in allocate_memory, or other parts in the code

	free(cluster.size_histogram.bin);
	free(cluster.size_distribution.bin);

	free(&copyslice[0]);
	free(&start_slice[0]);
	free(&slice[0]);

	if (sys.sim_type==MC_ALGORITHM){
        free(&cp_slice[0]);
    }

	return;
}

void special_init(Slice *psl){
    // printf("set the reduced bond number matrix\n\n ");
    total_energy(psl);
    if (psl->nbonds==psl->nparts){        set_reduced_bondnumber_ring(psl);}
    else if(psl->nbonds==psl->nparts-1){  set_reduced_bondnumber_chain(psl);}
    else{
        dprint(psl->nbonds);
        dprint(psl->nparts);
        error("something wrong with initializing the reduced bond numbers. the input must be a chain or ring structure");
    }

    // TST; hard coded lambda interface values...
    if (langevin.bond_breakage_analysis==1){
        // see Table 1  of  DOI: 10.1039/d3sm01255g
        tst_lambda.S_lambda_12=0.7; 
        tst_lambda.lambda_12=tst_lambda.S_lambda_12*pot.E_smin;

        tst_lambda.S_lambda_23=0.01;
        tst_lambda.lambda_23=tst_lambda.S_lambda_23*pot.E_smin;
        
        tst_lambda.bondbreakage_treshold=0.5;
    }

}