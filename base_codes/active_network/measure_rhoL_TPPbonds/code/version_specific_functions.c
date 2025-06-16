#include "path.h"
#include "analysis.h"


// ***************  LOCAL DENSITY ***************

/*------------------LOCAL FUNCTIONS------------------------------------------*/
void local_density(Slice *);
void N_TPP_bonds(Slice *);
/*---------------------------------------------------------------------------*/


void version_specific_analysis(Slice *psl){
	// this funciton is located inside terminate_block() at the outerloop in main.c	

    local_density(psl);
    N_TPP_bonds(psl);

}

void innerloop_analysis(Slice *psl){
	// this funciton is located inside in innerloop in main.c
	// so after ncycle2 bmd, mc or trajectory reads, you perform this function. 
}



void local_density(Slice *psl){
    int j = 0;
    double max_ld = 2.0, binwidth_ld = 0.038;
    int number_bins_ld = (int)ceil(max_ld / binwidth_ld);
    static int init = 1;

    for (int nbins = 2; nbins < 12; nbins++) {
        if (init) {
            if (j >= MAXBIN) {
                fprintf(stderr, "ERROR: j (%d) exceeds MAXBIN (%d)\n", j, MAXBIN);
                exit(EXIT_FAILURE);
            }

            ld_array[j].box_gridarea = (sys.boxl.x * sys.boxl.y) / (nbins * nbins);
            ld_array[j].ld_gridsize = max_ld / number_bins_ld;
            ld_array[j].local_density.length = number_bins_ld;

            ld_array[j].local_density.bin = (Statistics *)calloc(number_bins_ld, sizeof(Statistics));
            if (ld_array[j].local_density.bin == NULL) {
                fprintf(stderr, "Error allocating memory for ld_array[%d].local_density.bin\n", j);
                exit(EXIT_FAILURE);
            }

            for (int i = 0; i < number_bins_ld; i++) {
                ld_array[j].local_density.bin[i].mean = 0.0;
                ld_array[j].local_density.bin[i].variance2 = 0.0;
                ld_array[j].local_density.bin[i].n = 0;
            }
        }

        int *grid_count = (int *)calloc(nbins * nbins, sizeof(int));
        int *ld = (int *)calloc(number_bins_ld, sizeof(int));
        if (!grid_count || !ld) {
            error( "Memory allocation failed for grid_count or ld\n");
        }

        double gridsize_x = sys.boxl.x / (double)nbins;
        double gridsize_y = sys.boxl.y / (double)nbins;

        for (int ipart = 0; ipart < psl->nparts; ipart++) {
            vector pos = psl->pts[ipart].r;
            pbc(pos, sys.boxl);

            if (fabs(pos.x - 0.5 * sys.boxl.x) < 1e-6) pos.x = -0.5 * sys.boxl.x;
            if (fabs(pos.y - 0.5 * sys.boxl.y) < 1e-6) pos.y = -0.5 * sys.boxl.y;

            double x = (pos.x + 0.5 * sys.boxl.x) / gridsize_x;
            double y = (pos.y + 0.5 * sys.boxl.y) / gridsize_y;

            int binx = (int)x;
            int biny = (int)y;
            int bin_grid = binx * nbins + biny;

            if (x < 0 || binx >= nbins) gprint(x);
            if (y < 0 || biny >= nbins) gprint(y);

            if (bin_grid >= nbins * nbins || bin_grid < 0) {
                dprint(ipart);
                gprint(psl->c_time);
                dprint(bin_grid); dprint(nbins);
                dprint(binx); gprint(x);
                gprint(psl->pts[ipart].r.x); gprint(pos.x); gprint(gridsize_x); gprint(sys.boxl.x);
                dprint(biny); gprint(y);
                gprint(psl->pts[ipart].r.y); gprint(pos.y); gprint(gridsize_y); gprint(sys.boxl.y);
                error("bin_grid larger than nbins*nbins or < 0");
            }

            grid_count[bin_grid] += 1;
        }

        double grid_area = gridsize_x * gridsize_y;
        int tot_count = 0;

        for (int a = 0; a < nbins * nbins; a++) {
            double count = (double)grid_count[a];
            tot_count += count;

            double local_density = count / grid_area;
            if (local_density > max_ld) continue;

            if (local_density > 2.0) {
                gprint(gridsize_x); gprint(sys.boxl.x); dprint(nbins);
                for (int ipart = 0; ipart < psl->nparts; ipart++) {
                    int binx = (int)((psl->pts[ipart].r.x + 0.5 * sys.boxl.x) / gridsize_x);
                    int biny = (int)((psl->pts[ipart].r.y + 0.5 * sys.boxl.y) / gridsize_y);
                    int bin_grid = binx * nbins + biny;
                    if (bin_grid == a) {
                        dprint(ipart);
                        vprint(psl->pts[ipart].r);
                    }
                }
                error("local density larger than 2. Probably because particles are stacking.");
            }

            int bin_local_density = (int)(local_density / ld_array[j].ld_gridsize);
            if (bin_local_density >= number_bins_ld || bin_local_density < 0) {
                gprint(local_density); gprint(ld_array[j].ld_gridsize); dprint(bin_local_density);
                error("bin_local_density >= number_bins_ld or < 0");
            }

            ld[bin_local_density] += 1;
        }

        if (tot_count != psl->nparts) {
            dprint(tot_count);
            error("totcount != nparts");
        }

        for (int l = 0; l < number_bins_ld; l++) {
            running_statistics(&ld_array[j].local_density.bin[l], ld[l]);
        }

        char filename[100];
        snprintf(filename, sizeof(filename), "local_density_A%5.5lf.out", grid_area);
        FILE *fp = fopen(filename, "w");
        if (!fp) {
            printf("Warning: cannot open %s\n", filename);
            error("stop");
        } else {
            for (int a = 0; a < number_bins_ld - 1; a++) {
                double phi1 = a * ld_array[j].ld_gridsize;
                double phi2 = (a + 1) * ld_array[j].ld_gridsize;
                fprintf(fp, "%.4lf %.4lf %.10lf %ld\n", phi1, phi2,
                        ld_array[j].local_density.bin[a].mean,
                        ld_array[j].local_density.bin[a].n);
            }
            fclose(fp);
        }

        free(grid_count);
        free(ld);
        j++;
    }

    init = 0;
    return;
}

void N_TPP_bonds(Slice *psl){
	// the number of bonds that the TPP particles have

	int i, npatches[2];
	npatches[0]=3;
	npatches[1]=2;

	//loop over particles and check if it is a TPP partcle
	for (i=0; i<2 ; i++){
		int nbonds_P=0,n_P=0,node_P=0,dibond_P=0,tribond_P=0, monomer_P=0, monobond_P=0;
		for (int ipart=0; ipart<psl->nparts;ipart++){
			// for TPP
			if (sys.particletype[psl->pts[ipart].ptype].nsites==npatches[i]){
				n_P+=1;
				nbonds_P+=psl->pts[ipart].nbonds;
				if (psl->pts[ipart].nbonds==3){
					tribond_P+=1;
				}
				else if (psl->pts[ipart].nbonds==2){
					dibond_P+=1;
				}
				else if (psl->pts[ipart].nbonds==1){
					monobond_P+=1;
				}
				else if (psl->pts[ipart].nbonds==0){
					monomer_P+=1;
				}
				else{
					error("not able to determine number of bonds TPP particle makes?");
				}
			}
		}
		double average_nbonds_P=(double)nbonds_P/n_P;

		char value[1000];
	
		snprintf(   value, sizeof( value ), "%.2lf %.3lf %d %d %d %d\n", psl->c_time,average_nbonds_P,tribond_P, dibond_P, monobond_P, monomer_P  );

		if (i==0) write_append_to_file("N_TPP_bonds",  "", 'a', value );
		else if (i==1) write_append_to_file("N_DP_bonds",  "", 'a', value );
	}

    return;
}


/*_________________________MEMORY ALLOCATION AND FREEING________________*/

void special_init(Slice *psl){
    // for cluster analysis 
    if(cluster.analysis==1){
        // chain. = (Statistics *)calloc(NPART,sizeof(Statistics));
        printf("\nThe code performs cluster analysis \n");
        printf("      the total energy is %lf\n", total_energy(psl));
        printf("      now cluster analysis\n");
        psl->nclusters = cluster_analysis(psl); 
        clustersize_identification(psl);

        printf("      Initial # bonds is %d \n", psl->nbonds);
     
        strcpy(cluster.size_histogram.filename,"clustersize_histogram.out");
        strcpy(cluster.size_distribution.filename,"clustersize_distribution.out");
 
    }
    
}
void allocate_memory(Slice *psl){
	// allocate versions specific memory

	return;
}

void free_all_memory(void){
	// free the memory you allocated in allocate_memory, or other parts in the code


	free(&start_slice[0]);
	free(&slice[0]);


	for (int j = 0; j < 10; j++) {
		free(ld_array[j].local_density.bin);
		ld_array[j].local_density.bin = NULL;
	}
	return;
}
