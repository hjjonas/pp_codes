#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "path.h"

// ############### all function special to active_network.c ############################
// calculating time correlation of dissociation and association. 
// you can use read_coordinates_from_file(Slice *current, int block_number) to read trajectories.
// use int block_number as tau 
// slide over the data with 10 second timegaps (tau) and run over all sites that have not terminated yet
// so you need to identify whether sites are free or bound first  (make a list)
// then loop over the bonds and free sites of the new loaded frame. 
// the bonds come in particle and site pairs
// the free sites are alone.

void setup_type2(Slice *);
int return_d_type2(Slice *, int ,  int  );
int return_reaction_type( int ,int  , int , int );

void perform_autocorrelation_measurement( int max_tau){
    int tau_0,tau_1,tau,site_n,bond_n,makesbond;
    int N_freesites,reaction_type;
    int isite,ipart,jpart,jsite,iptype,jptype;
    int itype2,jtype2; // second dimension of auto_dissociation_tau0 etc
    Slice *slice_tau0,*slice_tau1;

    Pts *psi1;

    // there are 3 types of binding/unbinding  reactions 
    // 1) DP + DP -> DP-DP
    // 2) TPP* + DP -> TPP*-DP   here the active force points TOWARD the bond formed
    // 3) TPP* + DP -> *TPP-DP   here the active force points AWAY   the bond formed
    // 4) all of the above

    for (tau_0=0;tau_0<max_tau;tau_0++){ // first slice
        // make the free site and bond list now
        // printf("make the free site and bond list now in &slices[%d]\n",tau_0);
        // gprint(slices[tau_0].energy);
		slice_tau0=&slices[tau_0];

        identify_bound_and_freesites(slice_tau0);
        setup_type2(slice_tau0);
        // printf("\n \n");
        printf("  >>>>>>>>>>>>Nfree=%zu Nbond=%zu in tau_0=%d <<<<<<<<\n", free_sites_ipart->used,bondinfoarray->used,tau_0);
        // printf("  auto_association_tau0[ptype][type2][tau] \n");

        

        for (tau = 0; tau < max_tau-tau_0; tau++) {  //increment tau
            tau_1= tau_0+tau;
            slice_tau1=&slices[tau_1];
            // loop over FREE SITES, are they still free? 

            // printf("\ntau_1=%d  (tau_1= tau_0+tau= %d + %d) \n", tau_1,tau_0,tau);
            // printf("loop over all free sites (%zu) \n",free_sites_ipart->used);
            // printIntsArray(free_sites_ipart);
            // printBondTimeInfoArray(bondinfoarray );

            // every reaction decay is  measures for each tau, first normalize to probability, then add to 
            for( site_n=0 ; site_n<free_sites_ipart->used; site_n++){
            	// walk over the free sites, that are indicated by an ipart and isite, 
            	// free_sites_ipart->used is lowered when it makes a bond, so in the end there is no free site to loop over

                isite=free_sites_ipart->isite[site_n];
                ipart=free_sites_ipart->ipart[site_n];
                
                // and a ptype (DP or TPP) and type2 (toward or away from active force for TPP, no difference for DP)
                iptype=slice_tau0->pts[ipart].ptype;      // should this be of the tau_0 
                itype2=slice_tau0->pts[ipart].a_type2[isite];    // it should be the type2 of tau0            
                
                makesbond=0;
                // check is the slice_tau1 if this patch makes a bond 
                psi1=&slice_tau1->pts[ipart]; 

                // printf("   ipart=%d, type%d type2%d isite=%d  with ipart Nbond=%d\n",ipart,ptype,type2,isite,slices[tau_1].pts[ipart].nbonds);
                // loop over the bonds of ipart and see if isite is there
                for (bond_n=0; bond_n<psi1->nbonds;bond_n++){
                    if (isite==psi1->bound_site[bond_n]) {
                        makesbond=1; // yes it makes a bond now

                        jpart= psi1->bonds[bond_n];
                        // printf("**** \n");
                        // printf("   ipart %d and jpart %d make a bond \n", ipart, jpart);
                        // printf("   jpart makes bonds with: \n");

                        for (int bondj=0; bondj<slice_tau1->pts[jpart].nbonds ; bondj++) {
                        	// dprint(slice_tau1->pts[jpart].bonds[bondj]);
                        	// printf("  now find jsite out of %d bonds :\n",slice_tau1->pts[jpart].nbonds);
                        	// dprint(slice_tau1->pts[jpart].bonds[bondj]);
                        	if (ipart==slice_tau1->pts[jpart].bonds[bondj]){
                        		// printf("  the last one it is. \n");
                        		
	                        	jsite=slice_tau1->pts[jpart].bound_site[bondj];
	                        	break;
	                        }
                        }

                        // printf("   !!!ipart %d isite %d creates a bond with jpart=%d jsite=%d at tau=%d!!\n",ipart,isite,jpart,jsite,tau);
                        break;
                    }
                }

                // if it makes a bond, then remove it from the free_sites_ipart list 
                // else , plus 1 for the autocorrelation of association at time tau
                if (makesbond){
                  
                    removeElementXIntsArray(free_sites_ipart, ipart, isite );
                    // lower sites_n, you just removed isite, and now the next site is in the same posistion in the list
                    site_n--;
                	// printf(" (ipart/isite %d / %d)  made bond, %zu free bonds left \n",ipart,isite,free_sites_ipart->used);
				    // printf(" the next particle is (ipart/isite %d / %d) at position site_n=%d \n",free_sites_ipart->ipart[site_n],free_sites_ipart->isite[site_n] ,site_n);
                   
                    // now identify what type of reaction happened (reaction 1, 2 or 3 see above)
                    // check if jpart/jsite is in the free_sites_ipart, if yes, remove it there too
                    // jpart/jsite could be available from a previously broken bond between [tau0,tau1]
                    if (checkElementXIntsArray(free_sites_ipart, jpart, jsite)) removeElementXIntsArray(free_sites_ipart, jpart, jsite );

	            	// ok, ipart and jpart make a bond. but which reaction type is it? 
	            	jptype=slice_tau0->pts[jpart].ptype;       
	        		jtype2=slice_tau0->pts[jpart].a_type2[jsite];        
					
					reaction_type= return_reaction_type(  iptype,  itype2,  jptype,  jtype2);

					// so I need 3 association histograms (for the 3 reactions). make them in main.c
					// dprint(reaction_type);

					// mmm now I only will count if an association happens. does not look at particles that do not associate at all. 
					for (int teller=0; teller<tau; teller ++){
						auto_association_tau0[reaction_type][teller]+=1;
						auto_association_tau0[3][teller]+=1; // always count total reactivity
					}
                } 
                else{
                    // printf("   auto_association_tau0[ptype][type2][tau]+=1 with tau=%d (still free patch) and ptype %d, type2 %d\n",tau,ptype,type2);
                    if( (iptype>2) || (itype2>5) || (itype2<0)|| (tau>max_tau)){
                		error("(iptype>sys.nparticle_types) || (itype2>2) || (tau>max_tau)");
                	}
                }
            }

            // loop over the BONDS
            // printf("loop over all bonds (%zu) \n",bondinfoarray->used);
            for( bond_n=0 ; bond_n<bondinfoarray->used;bond_n++){
            	// bonds are made up by two particles and their corresponding patches 
                ipart=bondinfoarray->bondinfos[bond_n].ipart;
                iptype=slice_tau0->pts[ipart].ptype;
                itype2=bondinfoarray->bondinfos[bond_n].itype2;

                jpart=bondinfoarray->bondinfos[bond_n].jpart;
                jptype=slice_tau0->pts[jpart].ptype;
                jtype2=bondinfoarray->bondinfos[bond_n].jtype2;

                if ((jtype2<0) || (itype2<0) ){
                	error("(jtype2<0) || (itype2<0)"); 
                }

                // are they broken in slice_tau1? i.e. r>0.5
                // with the same distance criterion as for breakage paper 
                if( particles_distance(slice_tau1,   ipart,  jpart)>=0.5){
                    // they broke; remove them from bondinfoarray
                    // printf("   !!!particle %d and %d broke remove them from bondinfoarray\n",ipart,jpart);
                    removeIthElementBondTimeInfoArray(bondinfoarray,bond_n);
                    bond_n--;

                    reaction_type=return_reaction_type(  iptype,  itype2,  jptype,  jtype2);

                    // printf("theres is a dissociation reactoin %d detected at tau=%d\n",reaction_type,tau );
                    for (int teller=0; teller<tau; teller ++){
						auto_dissociation_tau0[reaction_type][teller]+=1;
						auto_dissociation_tau0[3][teller]+=1; // always count total reactivity

					}
					
                }
                else{
                    // still bound, nothing changes so +1
                 	if( (iptype>sys.nparticle_types) || (itype2>5) || (tau>max_tau)){
                		error("(iptype>sys.nparticle_types) || (itype2>2) || (tau>max_tau)");
                	}
                }
            }

        }

        //  it there are free patches or particles left, add then to the auto_tau0 
        for( int t=0; t<max_tau-tau_0; t++){
        	auto_association_tau0[3][t]+=(int)free_sites_ipart->used;
        	auto_dissociation_tau0[3][t]+=(int)bondinfoarray->used;
        }

        // is this necessary? I wont get std devs but im only interested in the average anyway.
        int reaction_types=4;
        for (int n=0; n<reaction_types;n++){
            update_autocorrelationfunction(&auto_association[n] , auto_association_tau0[n] , max_tau-tau_0);
            update_autocorrelationfunction(&auto_dissociation[n] ,auto_dissociation_tau0[n] ,max_tau-tau_0);	        
        }
    }
    // error("stop");


    return;    
}


int return_reaction_type( int itype,int  itype2, int jtype, int jtype2){
	int check_type, check_type2;
	// first check which particle is the DP particle, then check the other particle
	if (itype==1){
		check_type=jtype;
		check_type2=jtype2;
	}
	else if (jtype==1){
		check_type=itype;
		check_type2=itype2;
	}
	else{
		// if none of the particles are a DP (always a second particle type)
		// maybe the second particle type is not DP
		dprint(itype);
		dprint(jtype);
		error("TPP (particletype 0) and TPP  (particletype 0) are making bond?");
	}

	// reaction 1)
    if (check_type==1  ) return 0; // other particle also DP 
    else if ( (check_type==0) && (check_type2==0)) return 1;
    else if ( (check_type==0) && (check_type2==1)) return 2;
    else{ error("not able to determine if addition reaction is 1, 2, or 3 ");
	}
	return -1;

}

void setup_type2(Slice *psl){
	Pts *psi;
	

    vector Fa;
    tensor rotmat;

	//association type2
	for (int ipart=0; ipart<psl->nparts;ipart++){
		psi=&psl->pts[ipart];

        // important for the association ; 
		if (sys.particletype[psi->ptype].nsites==2){ //DP
			psi->a_type2[0]=0; psi->a_type2[1]= 0;
			// if (psi->nbonds==0){ 		psi->a_type2[0]=0; psi->a_type2[1]= 0;}  
	        // else if (psi->nbonds==1){	psi->a_type2[0]=1; psi->a_type2[1]= 1;}  
	        // else if (psi->nbonds==2){	psi->a_type2[0]=-1; psi->a_type2[1]= -1;}  
	        // else {						psi->a_type2[0]=-1; psi->a_type2[1]= -1;}
	     // shouldnt be possible to look at dissociation if you are a monomer, gives error in perform_autocorrelation_measurement
	    }
	    else if (sys.particletype[psi->ptype].nsites==3){ // TPP
	    	// look at first patch vector
	    
	    	if  (sys.particletype[psi->ptype].activity==0){
	    		for (int isite=0; isite<sys.particletype[psi->ptype].nsites;isite++){
	    			psi->a_type2[isite]= 0;
	    		}
	    	}
	    	else{
	    		tensor rotmat = getrotmatrix(psi->q);
				matrix_x_vector(rotmat,sys.particletype[psi->ptype].active_force.r,Fa); 
				
				// this should be per patch, it depends on the bond.. 
				// vprint(Fa);
				// vprint(psi->patchvector[isite]);
				// gprint(vector_inp(Fa, psi->patchvector[isite]));

				for (int isite=0; isite<sys.particletype[psi->ptype].nsites;isite++){
			    	if (vector_inp(Fa, psi->patchvector[isite])>0){ // toward patch isitse
			    		psi->a_type2[isite]= 0;
				    	// if (psi->nbonds==0){ 		psi->a_type2[isite]= 0;}  
				        // else if (psi->nbonds==1){	psi->a_type2[isite]= 1;}  
				        // else if (psi->nbonds==2){	psi->a_type2[isite]= 2;} 
				        // else if (psi->nbonds==3){	psi->a_type2[isite]= -1;} 
			    	}else{		// away from patch isitse
			    		psi->a_type2[isite]= 1;
				    	// if (psi->nbonds==0){ 		psi->a_type2[isite]= 3;}  
				        // else if (psi->nbonds==1){	psi->a_type2[isite]= 4;}  
				        // else if (psi->nbonds==2){	psi->a_type2[isite]= 5;} 
				        // else if (psi->nbonds==3){	psi->a_type2[isite]= -1;} 
			        }
			    }
			}
	    }
	}

	// important for the dissociation. you (a patch) can only dissociate if you are part of a bond 
	for( int bond_n=0 ; bond_n<bondinfoarray->used;bond_n++){
		bondinfoarray->bondinfos[bond_n].itype2=return_d_type2(psl,  bondinfoarray->bondinfos[bond_n].ipart,   bondinfoarray->bondinfos[bond_n].isite );
		bondinfoarray->bondinfos[bond_n].jtype2=return_d_type2(psl,  bondinfoarray->bondinfos[bond_n].jpart,   bondinfoarray->bondinfos[bond_n].jsite );
	}
	return;
}

int return_d_type2(Slice *psl, int ipart,  int isite ){
 	//dissociation type2
	Pts *psi;
	psi=&psl->pts[ipart];

    vector Fa;
    tensor rotmat;

	if (sys.particletype[psi->ptype].nsites==2){ //DP
		return 0;
		// if (psi->nbonds==1){ 		return 0;}  
        // else if (psi->nbonds==2){	return 1;}  
        // // else if (psi->nbonds==3){	return 2;}  
        // else {						return -1;}
     // shouldnt be possible to look at dissociation if you are a monomer, gives error in perform_autocorrelation_measurement
    }
    if (sys.particletype[psi->ptype].nsites==3){ //TPP
    	// look at first patch vector
    	if  (sys.particletype[psi->ptype].activity==0){
    		return 0;
    	}
    	tensor rotmat = getrotmatrix(psi->q); 
		matrix_x_vector(rotmat,sys.particletype[psi->ptype].active_force.r,Fa); 
		// vprint(Fa);
		// vprint(psi->patchvector[isite]);
		// gprint(vector_inp(Fa, psi->patchvector[isite]))
		// this should be per patch, it depends on the bond.. 
    	if (vector_inp(Fa, psi->patchvector[isite])>0){ // beteen two patches
    		return 0;
	    	// if (psi->nbonds==1){ 		return 0;}  
	        // else if (psi->nbonds==2){	return 1;}  
	        // else if (psi->nbonds==3){	return 2;} 
    	}else{
    		return 1;		
	    	// if (psi->nbonds==1){ 		return 3;}  
	        // else if (psi->nbonds==2){	return 4;}  
	        // else if (psi->nbonds==3){	return 5;} 
        } // F_A toward  patch
    }
    return -1;
}

void update_autocorrelationfunction(StatsLength *a, unsigned int *b, int max_tau){

    // first normlaize b and add it to a 
    // printf("the first value of b is %d and max_tau=%d\n",b[0],max_tau);
    if (b[0]==0){
    	for(int l=0;l<a->length;l++){
	        b[l]=0; // reset b to zeros
	    }
    	return; // there were no free patches or bonds for this particle type. Happens for Fa=250 and -1 0 0 force direction.
    }

    for(int l=0;l<max_tau-1;l++){	 
        if (b[l]<b[l+1]){
        	dprint(l);
        	dprint(b[l]);
        	dprint(b[l+1]);
        	error(" probablity is increasing??  b[l]<b[l+1]");
        }
    }

    for(int l=0;l<max_tau;l++){	 
        running_statistics(&a->bin[l], (double) b[l]/ (double) b[0]);
    }
    for(int l=0;l<a->length;l++){
        b[l]=0; // reset b's to zeros
    }
   	return;

}

void identify_bound_and_freesites(Slice *current){
	// perform an energy calculation before this functions, to make the bond information available
	// only do this for the first slice at t=0
	// loop over all particles and its sites
	// 
	int isite, n, isite2,jpart,jsite,bound;
	BondTimeInfo bondinfo;  // a single bondinfo
 	
 	static int init=1;

 	if (init){
 		//initialize the free site and isite integer lists; dont forget to free them
 		printf("initialize the free site and isite integer lists. there are %d particles and %d bonds\n",current->nparts ,current->nbonds );
 		
 		//initialize the bond info array; dont forget to free them
 		bondinfoarray = (BondTimeInfoArray *)calloc(1,sizeof(BondTimeInfoArray));
 	    initBondTimeInfoArray( bondinfoarray,  current->nbonds+100) ; // +100 little margin

 	    // printf("initIntsArray( free_sites_isite,current->nparts); \n");
 	    free_sites_ipart = (IntsArray *)calloc(1,sizeof(IntsArray));

 		initIntsArray( free_sites_ipart, current->nparts+100); 

 	    init=0;
 	}

        // start all over with identifying bonds and free sites;
    bondinfoarray->used=0;
    free_sites_ipart->used=0;

	// printf(" used=%zu\n",free_sites_ipart->used);
	// print_slice_information(current);

	int teller=0;

 	// printf("\n>>>>  teller   IPART(type)   ISITE     bound/free?<<<<<\n");
	for (int ipart=0;ipart<current->nparts; ipart ++){
		//loop over sites
		for (isite=0; isite<sys.particletype[current->pts[ipart].ptype].nsites;isite++){

			// printf("    %d      %d (%d)      %d/%d    ",teller, ipart,current->pts[ipart].ptype, isite , sys.particletype[current->pts[ipart].ptype].nsites );
			bound=0; // not bound (yet)
			// if ipart makes bonds, check if isite makes bond, if not, add it to free sites
			if (current->pts[ipart].nbonds>0){
				for (n=0; n<current->pts[ipart].nbonds;n++){
					if( current->pts[ipart].bound_site[n]==isite){
						// this site is bound to jpart=psl->pts[ipart].bonds[n]
						bound=1;
						break;
					}
				}
			}
			if (bound){ // if isite is bound
				jpart=current->pts[ipart].bonds[n];
			
				if (ipart<jpart){ // take unique bonds, thus ipart<jpart
					// printf("    bound\n" );
					jsite=-1;
					for (n=0; n<current->pts[jpart].nbonds;n++){ // search for the jsite
						if( ipart==current->pts[jpart].bonds[n]){
							jsite=current->pts[jpart].bound_site[n];
							break;
						} 
					}
					if (jsite==-1) error("jsite=-1 ??");

					bondinfo.start_time=current->c_time;
	                bondinfo.end_time=0;
	                bondinfo.ipart=ipart;
	                bondinfo.isite=isite;
	                bondinfo.jpart=jpart;
	                bondinfo.jsite=jsite;

					insertBondTimeInfoArray(bondinfoarray, bondinfo);
				}
				// else{
				// 	// printf("    bound (already tracked)\n" );
				// }
			}
			else{
				insertIntsArray(free_sites_ipart, ipart, isite);
				// insertIntsArray(free_sites_isite, isite);
			}
			teller++;
		}
	}

	// check if the number identified of free patches is equal to the total number of patches 
	int maxsites=0;
	for (int n=0;n<sys.nparticle_types;n++){
		maxsites+=(sys.particletype[n].nsites*sys.particletype[n].nparticles);
	}
	if ((free_sites_ipart->used>maxsites) || (free_sites_ipart->used + bondinfoarray->used *2 != maxsites)){
		printf(" free_sites_ipart->used=%zu, bondinfoarray->used=%zu, maxsites=%d\n",free_sites_ipart->used,bondinfoarray->used, maxsites);
		error("free_sites_ipart->used>maxsites");
	}
	return;
}


void print_autocorrelations_tofile(void){
	int reaction_types=4;
	

	for (int n=0; n<reaction_types;n++){
		
		snprintf( auto_association[n].filename, sizeof( auto_association[n].filename ),  "association_reaction%d",n);
		snprintf( auto_dissociation[n].filename, sizeof( auto_dissociation[n].filename ),"dissociation_reaction%d",n);
	
	
		print_statistics_file(&auto_association[n], max_tau);
		print_statistics_file(&auto_dissociation[n], max_tau);
			
	}


	return;
}


// read the trajectory
int read_coordinates_from_file(Slice *current, int block_number){
    // the file vi bondbreakage_coordines.dat  is read
    // each block, separated by a white line, contains the folloing info:
    // 1: time, 2: particle number, 3-5: position, 6-9: quaternion
    // which are of course NPART lines per block

    static FILE *fp;
    static int measure_dtime=1,first=1,n_shots=0;
    static double dtime=0;

    int bufferLength=200;
    char str[bufferLength];
    int ipart,read_lines=0,WHITELINE=0;
    static int linecount=0,linecount_tracker=0, init=1;

    double time, x,y,z,time0;
    vector r;
    quaternion q;

    char filename[600];
	snprintf( filename, sizeof( filename ), "%s/trajectory.xyz", sys.directorypath);


    // printf("read_coordinates_from_file .. block_number=%d \n",block_number);
    if (init){
    	// open only once in the beginning
	    printf("read_coordinates_from_file %s\n",filename);
    	fp = fopen(filename,"r");
    	init=0;	
    }
    
    
    if(fp == NULL) {
      error("Error opening file: trajectory.xyz");
      return -1;
    }
    // dprint(current->nparts+WHITELINE);
    current->nparts=slice[0].nparts;
    while(fgets(str, bufferLength, fp)!=NULL ) {
        if ((int)(linecount/(current->nparts+WHITELINE)) == block_number){
            //read the block and save it to the slice
            // printf(" reading snapshot number %d \n",block_number);
            // printf("the string: %s",str);
	        // dprint(linecount);
	        // dprint((int)(linecount/(sys.npart )));
	        // dprint(block_number);

            sscanf(str,"%lf %d %lf %lf %lf %lf %lf %lf %lf\n",&time,&ipart,&r.x,&r.y,&r.z,&q.q0,&q.q1,&q.q2,&q.q3);
            
            current->pts[ipart].r=r;
            current->pts[ipart].q=normalize_quat(q);
            // printf("read: %.1lf %d %lf %lf %lf %lf %lf %lf %lf\n",time,ipart,r.x,r.y,r.z,q.q0,q.q1,q.q2,q.q3);


            //update the patch vectors;
            current->pts[ipart].ptype=start_slice->pts[ipart].ptype; // WAS BUG!  Set manually that particle ipart has the same ptype as ipart from the start slice
            update_patch_vector_ipart(current,ipart); 
            
            if (read_lines==0){
            	// if (measure_dtime){
	            // 	if (first){  time0=time;	   first=0; }
	            // 	else{        dtime=time-time0; measure_dtime=0;}
	            // } 
	            // the time between two frames is always 1 second !! hardcoded
                current->c_time=block_number*1.0;
                // gprint(time);
                // gprint(current->c_time);
            }

            //check for errors
            if (read_lines!=ipart){
                printf("the %dth particle that is read is not ipart=  %d at time=%.1lf \n",read_lines,ipart,current->c_time);
                fclose(fp); // close the file
                error("ERROR: there is a problem in read_coordinates_from_file.");     
            }
            
            read_lines++;           
        }
        linecount++;
        if (read_lines==sys.npart){
            // printf("done reading block_number %d \n\n",block_number);
            // reset_center(current);
            // fclose(fp);
            n_shots+=1; // ?? not used
            // dprint(n_shots);
            return 1;
        }
    }

    if (fgets(str, bufferLength, fp)==NULL){
    	fclose(fp); // close the file
    	// error(" end of file of trajectory.");
    	return -1;
    }
    
    return 1;

}



void local_density(Slice *psl){
	/*   Local Density
     When MIPS (Motility Induced Phase Separation) occurs, there are high and low density regions.

     make a grid , and count number of particle in that grid --> local density
     for N=500, the dimenstions are  36.18 x 36.18 
     misschien 10x10 grid and/or 15x15

     de local density is a 1 one-dimensional array, that is a global variable ld_array
	  */


	int j=0;
	double max_ld=2.0,binwidth_ld=0.038;
	int number_bins_ld=(int)ceil(max_ld/binwidth_ld);

	static int init=1;

	//select a few number if grid bins, because I dont know how many I need
	for (int nbins=2 ; nbins<12; nbins++){ // nbins = 2,3,4,...,8 
		if (init){  
			ld_array[j].box_gridarea=(sys.boxl.x*sys.boxl.y)/(nbins*nbins); //the area of the grid that divides the box
			ld_array[j].ld_gridsize=(double)(max_ld/number_bins_ld);
		}  // hardcoded to 50 bins and range [0,2]  // I'm not sure if the ld_array[j].local_density array contains zeros

		// grid_count as 1D
		int *grid_count;
		grid_count =(int *) calloc(nbins*nbins, sizeof(int));
		
		for (int i=0 ; i<nbins*nbins ; i++){
			grid_count[i]=0;
		}

		int *ld;
		ld =(int *) calloc(number_bins_ld, sizeof(int));
		// for (int l=0; l<number_bins_ld; l++){
		// 	dprint(ld[l]);
		// }
		double gridsize_x=sys.boxl.x/(double)nbins;
		double gridsize_y=sys.boxl.y/(double)nbins;
		
		// loop over particles, and put them in a bin
		for (int ipart=0;ipart<psl->nparts;ipart++){
			vector pos=psl->pts[ipart].r;

			pbc(pos,sys.boxl);
			// ok. due to the printing accuracy of the trajectory file. it may happen (quite often) that the x or y position is exactly on th boundary
			//so pbc won't catch it because a.x==b.x. (only positive values of x and y, negative values will lead to bin =0 )
			
			if ( fabs(pos.x - 0.5*sys.boxl.x )< 1e-6) { // i.e. if pox.x==0.5*box.x
				pos.x=- 0.5*sys.boxl.x ;
			}
			if ( fabs(pos.y - 0.5*sys.boxl.y )< 1e-6) { // i.e. if poy.y==0.5*boy.y
				pos.y=- 0.5*sys.boxl.y;
			}

			double x=(pos.x+0.5*sys.boxl.x)/gridsize_x ;
			double y=(pos.y+0.5*sys.boxl.y)/gridsize_y ;
			int binx=(int)x ;
			int biny=(int)y  ;
			
			int bin_grid=binx*nbins + biny;
			if (x<0 ||  binx>=nbins){ gprint(x);
		    }
		    if (y<0||  biny>=nbins){ gprint(y);
		    }
			if (bin_grid>nbins*nbins || bin_grid<0){ 
				dprint(ipart);
				gprint(psl->c_time);
				dprint(bin_grid);   dprint(nbins);
				dprint(binx);
				gprint(x);
				gprint(psl->pts[ipart].r.x);
				gprint(pos.x);
				gprint(gridsize_x);
				gprint(sys.boxl.x);
				printf(" \n");

				dprint(biny);
				gprint(y);
				gprint(psl->pts[ipart].r.y);
				gprint(pos.y);
				gprint(gridsize_y);
				gprint(sys.boxl.y);

				error("bin_grid largen than nbins*nbins or < 0");}

			grid_count[bin_grid ]+=1; // add one to the grid count 1D
			
		}
		

		// loop over grid_count, calculate the local density in the grid and add it to 
		// local density is max 2 , the binsize of local density is 2/100
		double grid_area=gridsize_x*gridsize_y;

		int tot_count=0; 
		for (int a=0;a<nbins*nbins;a++){ //1D
			double count=(double)grid_count[a]; 
			tot_count+=count;

			//calculate local_denisty in the grid
			double local_density=count/grid_area;
			// double local_density=count/sys.npart*5.;

			if( local_density>max_ld){
				// printf("WARNING: local_density=%lf > max_ld=%lf",local_density,max_ld);
				continue;
			}
			if( local_density>2){  // error
				gprint(gridsize_x);
				gprint(sys.boxl.x);
				dprint(nbins);
				for (int ipart=0;ipart<psl->nparts;ipart++){

					int binx=(int)((psl->pts[ipart].r.x+0.5*sys.boxl.x)/gridsize_x)  ;
					int biny=(int)((psl->pts[ipart].r.y+0.5*sys.boxl.y)/gridsize_y ) ;

					
					int bin_grid=binx*nbins+biny;
					
					if (bin_grid==a){
						dprint(ipart);
						vprint(psl->pts[ipart].r);					
					}				
				}
				error("local density larger than 2. Probably because particles are stacking."); 
			} // error
			
			// add one to global array 
			int bin_local_density=(int)(local_density/ld_array[j].ld_gridsize);
			// int bin_local_density=(int)(local_density*number_bins_ld);

			if (bin_local_density>number_bins_ld || bin_local_density<0){
			  gprint(local_density);
			  gprint(ld_array[j].ld_gridsize);
			  dprint(bin_local_density); error("bin_local_density larger than number_bins_ld or < 0");
			}
			ld[bin_local_density]+=1; //just add one to the histogram

		}
		
		if (tot_count!=psl->nparts){
			dprint(tot_count);
			error("totcount != nparts");
		}
		//update statistics
		for (int l=0;l<number_bins_ld;l++){ //1D
			running_statistics(&ld_array[j].local_density.bin[l], ld[l]);
		}


		// print to file
		FILE *fp;
		char filename[100];
		snprintf( filename, sizeof( filename ), "local_density_A%5.5lf.out", grid_area);
	    if((fp=fopen(filename,"w"))==NULL) { // overwrite the file
	        printf("Warning: can not open %s\n",filename);
	        error("stop");
	    }
	    else {
	        for (int a=0;a<number_bins_ld-1;a++){ 
	        	double phi1=a *ld_array[j].ld_gridsize;
	        	double phi2=(a+1) *ld_array[j].ld_gridsize;
	  
	        	fprintf(fp,"%.4lf %.4lf %.10lf %ld\n", phi1,phi2,ld_array[j].local_density.bin[a].mean, ld_array[j].local_density.bin[a].n );
	        }
	    }
	    fclose(fp); //close the file

		//free the memory of the grid_count and ld
		free(grid_count); 
		free(ld); 

		j++;
	}
	init=0;
    return;
}

void N_TPP_bonds(Slice *psl){
	// the number of bonds that the TPP particles have

	int nbonds=0,ntpp=0,node=0,dibond=0,tribond=0, monomer=0, monobond=0;
	//loop over particles and check if it is a TPP partcle
	for (int ipart=0; ipart<psl->nparts;ipart++){
		if (sys.particletype[psl->pts[ipart].ptype].nsites==3){
			ntpp+=1;
			nbonds+=psl->pts[ipart].nbonds;
			if (psl->pts[ipart].nbonds==3){
				tribond+=1;
			}
			else if (psl->pts[ipart].nbonds==2){
				dibond+=1;
			}
			else if (psl->pts[ipart].nbonds==1){
				monobond+=1;
			}
			else if (psl->pts[ipart].nbonds==0){
				monomer+=1;
			}
		}
	}
	double average_nbonds=(double)nbonds/ntpp;

	// print to file
	FILE *fp;
	char filename[100];
	snprintf( filename, sizeof( filename ), "N_TPP_bonds.out");
    if((fp=fopen(filename,"a"))==NULL) { // overwrite the file
        printf("Warning: can not open N_TPP_bonds.out\n");
    }
    else {
        fprintf(fp,"%.2lf %.3lf %d %d %d %d\n", psl->c_time,average_nbonds,tribond, dibond, monobond, monomer );    
    }
    fclose(fp); //close the file

    return;
}

void print_bondtracking(BondTimeInfoArray *a,  int elemi , int broken_bound, double ctime){
    //broken_bound is an integer, its 0 if the bond is broken , or 1 if it has formed
    // if
    FILE *fp;
    double lifetime;

    if((fp=fopen("bondtracking.out","a"))==NULL) {; // you cannot call it bondtracking.out, because it might be deleted
        printf("Warning: can not open bondtracking.out\n");
    }
    else {
        if (broken_bound==-1){
        	//if the bond is broken, you know at which elemi the info is stored
        	// after printing to the files, you remove this info.
            lifetime=a->bondinfos[elemi].end_time-a->bondinfos[elemi].start_time;
        }else if (broken_bound==1){
            elemi=a->used-1; // why -1; because before  print_bondtracking, I added the new bond. So the last added bond is accessable via a->used-1
            lifetime=0;
        }
        fprintf(fp,"%d %3.2f %d %d %.3f\n", broken_bound, ctime,a->bondinfos[elemi].ipart,a->bondinfos[elemi].jpart,lifetime);
    }
        
    fclose(fp);
    return;
}

