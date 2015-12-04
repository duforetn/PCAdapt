/*
 *    FastPCAdapt Cov_line.c
 *    Copyright (C) 2014 Nicolas Duforet Frebourg
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    this program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "math.h"
#include "Cov_line.h"
#include "../src/matrix.h"
#include "../src/Data.h"
#include "Data__f.h"
#include "../src/linAlgebra.h"

#define PRECISION_FLOAT 16
#define NA 9
#define TOL .0001

//get_row
////read a block of lines (snp), scale it and store it

int get_row(double *Geno, FILE *GenoFile, int nIND, double *mean, double *SNPSd, int *pairwiseObs, int sc, int blocksize){

	float value;
	int snp, ind, na_tot = 0, na, i, j;
	double var;
	int miss[nIND];
	double *presentSNP = calloc(nIND*blocksize, sizeof(double));
	double *scratchObs = calloc(nIND*nIND, sizeof(double));
	//int *tpresentSNP = calloc(nIND*blocksize, sizeof(int));
	for(snp=0; snp<blocksize; snp++){
		ind = 0;
		var = 0;
		na = 0;
		*mean = 0;
		while(ind < nIND){
				if(fscanf(GenoFile, "%g", &value) != EOF){
        		                Geno[snp*nIND + ind] = (double) value;
					if(value != NA){
						*mean += (double) value;
						miss[ind] = 0;
						presentSNP[snp*nIND + ind] = 1;
	//					tpresentSNP[ind*blocksize + snp] = 1;
					} else { na++;}
				}
				ind++;
		}
		//Update the pairwise missing count.
/*		for (i=0; i<nIND; i++){
                        for (j=i; j<nIND; j++){
                                if (!(miss[i] || miss[j])) {pairwiseObs[i*nIND + j]++; if (i != j) pairwiseObs[j*nIND + i]++;}
                        }
                }
*/		if (nIND <= na)	{
			*mean = NA;
		} else {
			*mean /= (nIND - na);
		}
		for (ind=0; ind<nIND; ind++) { if(Geno[snp*nIND + ind] != NA) {Geno[snp*nIND + ind] -= *mean; var += Geno[snp*nIND + ind]*Geno[snp*nIND + ind];} else { Geno[snp*nIND + ind] = 0; } } 
		if (sc){
			if (var > TOL){
				var = var/(nIND - na);
				for (ind=0; ind<nIND; ind++) {Geno[snp*nIND + ind] /= sqrt(var);}
//mAF filtering
//				if (*mean < .1 || *mean > 1.9) for (ind=0; ind<nIND; ind++) {Geno[snp*nIND + ind] = 0;}
				SNPSd[snp] = sqrt(var);
			} else {
				SNPSd[snp] = 1;
			}
		} else {
			SNPSd[snp] = sqrt(var/(nIND - na));
		}
		na_tot += na;
	}	
	add_to_obs(pairwiseObs, scratchObs, nIND, presentSNP, blocksize);
	free(scratchObs);
	free(presentSNP);
	return na_tot;
}

//add_to_cov
//Given a block of snps, computes correlation matrix of the block and add it to the existing matrix.

void add_to_cov(double *Cov, double *scratchCov, int nIND, double *Geno, int blocksize){

	tAA(Geno, scratchCov, blocksize, nIND);
	int i;
	for (i=0; i<nIND*nIND; i++) Cov[i] += scratchCov[i];

}

void add_to_obs(int *pairwiseObs, double *scratchCov, int nIND, double *Geno, int blocksize){

        tAA(Geno, scratchCov, blocksize, nIND);
        int i;
        for (i=0; i<nIND*nIND; i++) pairwiseObs[i] += (int) scratchCov[i];

}


//Cov_line
//function to learn covariance (or cor) matrix, while reading a data file, by storing blocks of lines.
//TODO: parallelize

int Cov_line(double *Cov, double *SNPSd, int nSNP, int *nSNP_file, int nIND, int *pairwiseObs, int sc, char **GenoFileName, int nfile){
	FILE *GenoFile;
	double var, mean;
	int i, j, na, na_tot = 0, file;
	int *missing;
	int blocksize = 120, snp_count = 0;
	double *Geno = calloc(nIND*blocksize, sizeof(double));
	double *scratchCov = calloc(nIND*nIND, sizeof(double));

	for (file=0; file<nfile; file++){
		blocksize = 120;

	        if((GenoFile = fopen(GenoFileName[file], "r")) == NULL){
        	        printf("Error, invalid input file\n");
			return 1;
	        }
		for (i=0; i<nSNP_file[file] ; i += blocksize){
			if (nSNP_file[file] - i < blocksize) blocksize = nSNP_file[file] - i;
			na = get_row(Geno, GenoFile, nIND, &mean, SNPSd + i + snp_count, pairwiseObs, sc, blocksize);
			add_to_cov(Cov, scratchCov, nIND, Geno, blocksize);
			na_tot += na;
		}
		snp_count += nSNP_file[file];
		fclose(GenoFile);
	}

	if (na_tot) printf("%i out of %i missing data ignored\n", na_tot, nSNP*nIND);
	
	//Correct for missing data in each pair of individuals
	for (i=0; i<nIND; i++){
		for (j=0; j<nIND; j++){
			if (pairwiseObs[i*nIND + j] > 0){ Cov[i*nIND + j] *= (double) nSNP/pairwiseObs[i*nIND + j];}
	//		printf("%i ", pairwiseObs[i*nIND + j]);
		}
	//	printf("\n");
	}
	free(Geno);
	return 0;
}

