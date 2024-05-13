// version_specific_functions.h

#ifndef vsf_H
#define vsf_H
/*------------------GLOBAL VARIABLES-----------------------------------------*/

#define NBINS_MECH 25

/*---------------------------------------------------------------------------*/
/*------------------STRUCTURE DEFINITIONS------------------------------------*/



typedef struct TST_type {
// Transition State Theory, for each (reduced) bond
// see ref DOI: 10.1039/d3sm01255g                                      
// used in TST ; this is per simulation; after breakage it is reset in breakage_i_tst and saved in breakage_all_tst
  

  int               loc,  // previous, i.e. one timestep, region: below lambda_12 (loc=1), between lambda_12-lambda_23(loc=2), or above lambda_23 (loc=3)?
                    N_13,   // how often did you go from region 1 to 3 (first stage breakage)
                    N_31,   // how often did you go from region 3 to 1 (rebinding)
                    N_24,   // how often did you go from region 1 to 3 (second stage breakage)
                    N_1,   // how many positive crossings of lambda_12
                    flag_from1, // you come from  region 1 
                    flag_from3; // you come from  region 3  (use as 0 (off) and 1 (on))

                    
  double            last_lambda_12_time, // used for flux calculation
                    // these two variables are for tracking the mechanism. see eq 14-16 of the ref
                    tau_sim,           // when broken, this is the timestamp of the last lambda_23 crossing
                    S_mech,            // this is the tracked value of S of the mechanism. its sort of a dummy variable

                    // NOTE: these S histograms of length NBINS_MECH (see eq 14-16) are only used in breakage_all_tst as the distribution is indep whether the bond will break in the future!
                    S_hist_lambda_23_all[NBINS_MECH], // Fig 8b this is the histogram / distribution of all crossings through the lambda_23 interface
                    n_1to3[NBINS_MECH], // eq (14)
                    n_2to4[NBINS_MECH], // histogram for the mechanism as funtion of S
                    n_3to1[NBINS_MECH]; // histogram for the rebinding per S value


// average breakage properties, after many runs.
// the individual values are printed to a file after each breakage.
  Statistics          P_lambda23_lambda12,
                      Psep,
                      Prebind,
                      invflux_12,         // the average time (inverse flux) between two positive crossings of the lambda_12 interface 
                      tau_12;            // within one breakage event, how long does it take to go from lambda_12 to lambda_23 (on average)

} TST;

typedef struct TST_Lambda_type {
// defines the phase space as lambda (reaction coordinate), in this case, a function of (r,S)
// for all bonds the same

  double              bondbreakage_treshold,
                      S_lambda_12,  // value of the S at which the lambda_12 interface lies 
                      S_lambda_23,  // value of the S at which the lambda_23 interface lies 
                      lambda_12,   // value of the E_lambda at which the lambda_12 interface lies 
                      lambda_23;   // value of the E_lambda at which the lambda_23 interface lies 

} TST_Lambda;

/*---------------------------------------------------------------------------*/


/*------------------GLOBAL STRUCTURES-----------------------------------------*/
// the length of NPART is because the maximum number of bonds is NPART (for a ring)
extern TST breakage_all_tst[NPART],breakage_i_tst[NPART];
extern int bondnumbers[NPART][NPART], reduced_bondnumbers[NPART][NPART], Nred_bonds;
extern TST_Lambda tst_lambda;
extern Statistics tau_sims[NPART];
/*---------------------------------------------------------------------------*/


/*------------------GLOBAL FUNCTIONS-----------------------------------------*/
extern void innerloop_analysis(Slice *);
extern void version_specific_analysis(Slice *);
extern void every_timestep_analysis(Slice *);
extern void special_init(Slice *);
/*---------------------------------------------------------------------------*/


#endif
