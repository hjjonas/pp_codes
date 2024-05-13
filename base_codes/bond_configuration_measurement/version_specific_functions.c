#include "path.h"

//################################################################################################################
//######################################### BOND CONFIGURATION FUNCTIONS ##############################################
//################################################################################################################


/*
These version specific function measure the distribution inside the bond of a dimer, decamer or ring
it makes use of a 2D table reduced_bondnumbers[i][j] to lookup the (reduced) bondnumber which is made in init.c

Each BMDstep the S, r, and E are measured and histogrammed in histogramming_bond_configurations(psl)

*/

/*------------------LOCAL FUNCTIONS------------------------------------------*/
void histogramming_bond_configurations(Slice *);
int return_binnumber(double , double , double , double , int );
void update_histogram(Histogram *, int , int  );
void print_histogram( char [200], Histogram  );

void make_histogram(Histogram *, double , double , double );
void set_reduced_bondnumber_chain(Slice *);
void set_reduced_bondnumber_ring(Slice *);
/*---------------------------------------------------------------------------*/

Histogram  histogram_S, histogram_E, histogram_r;
int reduced_bondnumbers[NPART][NPART]={-2}, Nred_bonds;


void version_specific_analysis(Slice *psl){

    print_histogram( "all_bond_distances",  histogram_r );
    print_histogram( "all_bond_energies",  histogram_E );
    print_histogram( "all_switch_values",  histogram_S );

    return;
}

void every_timestep_analysis(Slice *psl){
    // perform these analysis function after every timestep. 
    // practically only used for bond breakage analysis or making histogrmas

    // see DOI: 10.1039/d3sm01255g
    // dprint(langevin.measure_bond_configurations);
    if (langevin.measure_bond_configurations){
        // printf("every_timestep_analysis\n");
        // eg Fig 8, 12, 16 in ref 
        // measures distributions when inside bond
        
        // first measure energy includes bond information
        psl->energy=total_energy(psl);

        if (psl->nbonds!=start_slice->nbonds ){
            // new configuration, and recalculate totalenergy
            restart_conf(psl);

            mc_warmup(psl);
            if (sys.nearest_neighbor){
                update_nnlist(psl);
            }
            psl->energy=total_energy(psl);
            if (psl->nbonds!=start_slice->nbonds ){ 
                dprint(psl->nparts);
                dprint(psl->nbonds);
                error("not enough bonds in conf.inp file.  shouldn't be possible");
            }
        }
        histogramming_bond_configurations(psl);
    }
}

void innerloop_analysis(Slice *psl){
	// this funciton is located inside in innerloop in main.c
	// so after ncycle2 bmd, mc or trajectory reads, you perform this function. 

    // If bond tracking is enabled, perform bond tracking for the current slice
    // active physical gels
    // if(analysis.bond_tracking==1){
    //     bond_tracking(psl);
    //     printBondTimeInfoArray(bondtimeinfoarray )
    // }

    return;
}


void histogramming_bond_configurations(Slice *psl){
    // loop over alle bodns, check which reduced bond it is and measure E, S, and r_ss. 
    // fill in the histogram based on these measurements
    // this code requires that bond infromation is updated. maybe add it to force, else perform energy calculation

    int redbond_nr,j,i,n,bin;
    double E,S,r_ss;

    // ####loop over all (unique) bonds ###
    for (i=0;i<psl->nparts; i++){
        for ( n=0;n<psl->pts[i].nbonds; n++){
            j=psl->pts[i].bonds[n];
            if (j<i){ continue ; } // // unique bodns
            redbond_nr=reduced_bondnumbers[i][j];
            // dprint(redbond_nr);
            // bond_energy
            bin=return_binnumber(psl->pts[i].bond_energy[n],histogram_E.min, histogram_E.max, histogram_E.d_bin, histogram_E.nbins);
            update_histogram(&histogram_E, bin, redbond_nr);
            
            // S value 
            bin=return_binnumber(psl->pts[i].bond_switch[n],histogram_S.min, histogram_S.max, histogram_S.d_bin, histogram_S.nbins);
            update_histogram(&histogram_S, bin, redbond_nr);
            
            //bond_distance 
            bin=return_binnumber(psl->pts[i].bond_distance[n],histogram_r.min, histogram_r.max, histogram_r.d_bin, histogram_r.nbins);
            update_histogram(&histogram_r, bin, redbond_nr);
        }
    }
   return;  
}

int return_binnumber(double value, double min, double max, double d_bin, int nbins){
    // is you fall outside the [min,max> range, you assing 0 or nbins-1, respectively. 
    //inside the range, the bin is determined by the d_bin, but as you start by 1, you add +1
    int bin_number;

    if (value<min ){         bin_number=0; }
    else if ( value>=max){   bin_number=nbins-1;}
    else{                    bin_number=(int)floor((value-min)/d_bin)+1;}
    
    if ((bin_number<0) || (bin_number>=nbins)){
        dprint(bin_number);
        dprint(nbins);
        error("bin_number out of range, <0 or >nbins");
    }
    return bin_number;
}

void update_histogram(Histogram *hist, int bin, int bondnr ){
   
    if ((bin>=0) && (bin<hist->nbins)){ // bin cannot be larger than hist->nbins
        // printf("update_histogram with count%d and bin %d",bondnr,bin);
        hist->bin[bondnr][bin]+=1;
        hist->length[bondnr]+=1;    
    }
    else{
        dprint(bin);
        dprint(hist->nbins);
        error("trying to update hist outside range");
    }
    return;
}

void print_histogram( char filename[200], Histogram hist ){

    int bufferLength=500;
    int length=0;
    FILE *fp;

    // print histograme to file: filename+"_histogram_"+length
    char filename2[bufferLength];
    // snprintf( filename2, sizeof( filename2 ), "%s_histogram_%ld",filename, hist.length[0]);
    snprintf( filename2, sizeof( filename2 ), "%s_histogram",filename);

    fp = fopen(filename2,"w");
    // printf("print histograme to file: filename+_histogram_+length\n");
    if (fp  == NULL){
        // printf("%s\n", &filename2);
        error("could not open file");
    }else{
        fprintf(fp, "#");
        for (int n=0; n<Nred_bonds; n++){
            fprintf(fp, "%ld ", hist.length[n]);
        }
        fprintf(fp, "\n ");
        // dprint(hist.nbins);
        for (int i=0; i<hist.nbins; i++){
            double lowerbound=hist.min + i*hist.d_bin;
            fprintf(fp, "%.5lf ", lowerbound);
            
            for (int n=0; n<Nred_bonds; n++){
                fprintf(fp, "%.15lf ", (double)(hist.bin[n][i])/hist.length[n]);
            }
            fprintf(fp, "\n ");
        }
    }
    fclose(fp);
    return;
}


/* *********** INITIALIZE STUFF**************/
 
 void special_init(Slice *psl){
    // printf("set the reduced bond number matrix ");
    total_energy(psl);
    if (psl->nbonds==psl->nparts){        set_reduced_bondnumber_ring(psl);}
    else if(psl->nbonds==psl->nparts-1){  set_reduced_bondnumber_chain(psl);}
    else{
        dprint(psl->nbonds);
        dprint(psl->nparts);
        error("something wrong with initializing the reduced bond numbers ");
    }

    make_histogram(&histogram_r,        0.0, 0.020,       0.020/50.);
    make_histogram(&histogram_E, pot.E_smin,   0.0, -pot.E_smin/50.);
    make_histogram(&histogram_S,        0.0,   1.0,           1/50.);

    printf("Setting up the histrograms for E,S,r is done\n");

    return;
}

void make_histogram(Histogram *hist, double min, double max, double d_bin){
    /*you need to make sure d_bin  yields an integer (of veeeeery close to it) */

   hist->min=(double)min;
   hist->max=(double)max;
   hist->d_bin=(double)d_bin;
   hist->nbins=(int)floor((hist->max-hist->min)/hist->d_bin)+2; 
   if ((hist->nbins>MAXBIN) || (hist->nbins==0)) error("nbins>MAXBIN  or ==0");

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
    for (int chain=0; chain<Nchains; chain++ ){
        // restart; set bond number to zero   
        int redbond_nr=0;
        // dprint(redbond_nr)
        // walk over chain of length Lchain
        for (int L=0; L<Lchain; L++){
            // printf(" the %d th particle = %d in chain %d of length %d\n",L,DPc,chain,Lchain );       
            // the bonds of DPc
            for (int n=0; n<psl->pts[DPc].nbonds;n++){

                int bound_particle=psl->pts[DPc].bonds[n];
                // printf("     particle %d and %d make bond  \n",DPc ,bound_particle);       

                // check if it as a DP particle int the list, or a TPP particle
                if ((checkElementXIntArray(&DP_list, bound_particle )==1) || (bound_particle<DP1) ){
                    // dprint(checkElementXIntArray(&DP_list, bound_particle ));
                    printf("      bond between (%d-%d) , given reduced bondnumebr:  %d \n",DPc,bound_particle,redbond_nr);       
                    // dprint(reduced_bondnumbers[DPc][bound_particle]);
                    reduced_bondnumbers[DPc][bound_particle]=redbond_nr;
                    reduced_bondnumbers[bound_particle][DPc]=redbond_nr;
                    // dprint(reduced_bondnumbers[DPc][bound_particle]);

                    if (bound_particle>= DP1){removeElementXIntArray( &DP_list ,  bound_particle);} 
                    redbond_nr+=1;
                }
            }
            if (checkElementXIntArray(&DP_list, DPc )==1){removeElementXIntArray( &DP_list ,  DPc);}
            DPc+=1;
        }
    }

    freeIntArray(&DP_list) ;
    // error("stop");
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
                red_bondnr++;
            }
        }
    }
    Nred_bonds=red_bondnr;
    if (Nred_bonds==0 ) error("no bonds found");
    dprint(Nred_bonds);
    return;
}