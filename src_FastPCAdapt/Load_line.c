/*
   FastPCAdapt Load_line.c
   Copyright (C) 2014 Nicolas Duforet Frebourg

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "../src/matrix.h"
#include "../src/Data.h"
#include "../src/linAlgebra.h"
#include "Data__f.h"
#include "math.h"
#include "Cov_line.h"
#include "Load_line.h"

#define PRECISION_FLOAT 16
#define NA 9

//Load_line
//computes the U matrix, by block of snps, with no storage.

int Load_line(double *U, double *Sigma, double *V, double *miss, double *mAF, int nSNP, int *nSNP_file, int K, int nIND, int *pairwiseObs, int sc, char **GenoFileName, int nfile){

	FILE *GenoFile;
	int i, ii, j, na, na_tot = 0, file, snp_count = 0;
	int *missing;
	double var, mean, ploidy = 2.0;
	double *Geno = calloc(nIND, sizeof(double));

        for (i=0; i<nIND; i++){
                for (j=0; j<K; j++){
                        V[i*K + j] /= Sigma[j];
                }
        }

	for(file=0; file<nfile; file++){
        	if((GenoFile = fopen(GenoFileName[file], "r")) == NULL){
	                printf("Error, invalid input file\n");
                	return 1;
        	}

		for (i=0; i<nSNP_file[file]; i++){
			ii = i + snp_count;
			na = get_row(Geno, GenoFile, nIND, &mean, &var, pairwiseObs, sc, 1);
			miss[ii] = na;
			mAF[ii] = mean;
			if(mAF[ii] != NA) mAF[ii] = (mAF[ii]/ploidy < 1 - mAF[ii]/ploidy) ? mAF[ii]/ploidy : 1 - mAF[ii]/ploidy;
			prodMatrix(Geno, V, (U + K*(ii)), 1, nIND, nIND, K);
		}

		snp_count += nSNP_file[file];
	        fclose(GenoFile);
	}

        for (i=0; i<nIND; i++){
                for (j=0; j<K; j++){
                        V[i*K + j] *= Sigma[j];
                }
        }

	free(Geno);
	return 0;
}

