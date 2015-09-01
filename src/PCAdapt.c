/*
   PCAdapt PCAdapt.c
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
#include <unistd.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "Data.h"
#include "matrix.h"
#include "../src_FastPCAdapt/FastPCAdapt.h"

int main(int argc, char* argv[]){

	if (argc == 1){
		Welcome(1);
		return 0;
	}
	if (!strcmp(argv[1], "fast")){
		return FastPCAdapt(argc, argv);
	}
	
        FILE *GenoFile;
	char *GenoFileName = calloc(256, sizeof(char));
	char *OutputFileName = calloc(256, sizeof(char));
	char *subSampleName = calloc(256, sizeof(char));
	int *missing;
	int nSNP, nIND, K = 2, runSVD = 0, tmp, prior = 1, sc = 1;
	int nsteps = 400, burnin = 200, na, transpose = 0;
	double alpha = .0001, Isingb = 0, prop = .01;
	double *Genotypes, *Factors, *Lambda, *Pz, *Piz, *sigma2_F;
	double Tau2 = 2, sigma2 = .1;
	double *sigma2_Lambda, *rowSd;
	double *C2;

	if (argc > 1){Welcome(0); if(handleParams(argc, argv, &nSNP, &nIND, &Genotypes, &GenoFileName, &OutputFileName, &K, &nsteps, &burnin, &runSVD, &subSampleName, &transpose, &prior, &Isingb, &sc, &prop)) return 0;} else {Welcome(1); return 0;}

	na = DataImputation(Genotypes, &missing, nSNP, nIND, "naive");
	if (na) printf("%i out of %i imputed data\n", na, nSNP*nIND);
	if(!strcmp(OutputFileName, "")) OutputFileName = "PCAdapt_output";
	printf("results in files %s*\n", OutputFileName);
	printf("Computation of %i factors\n", K);
	initializeVariables(&Factors, &sigma2_F, &Lambda, &Pz, &Piz, &C2, &sigma2_Lambda, &rowSd, K, alpha, nSNP, nIND);
	
	if (sc) scale(Genotypes, rowSd, nSNP, nIND, 1, sc);

	runMCMC(Genotypes, missing, Factors, sigma2_F, Lambda, Pz, Piz, &sigma2, K, alpha, Tau2, C2, sigma2_Lambda, Isingb, nSNP, nIND, nsteps, burnin, na, prior, runSVD, OutputFileName, subSampleName, prop);

	free(Pz);
	if(K > 1){free(Factors); free(Lambda);}
	free(sigma2_F);
	free(Piz);
	free(GenoFileName);
	free(Genotypes);
	free(missing);
	return 0;

}
