/*
   FastPCAdapt FastPCAdapt.c
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
#include "../src/Data.h"
#include "../src/matrix.h"
#include "../src/linAlgebra.h"
#include "Data__f.h"
#include "stat.h"

#define MAXFILE 64

int FastPCAdapt(int argc, char* argv[]){

        FILE *GenoFile;
	char **GenoFileName = malloc(sizeof(char)*MAXFILE);
	char *OutputFileName = malloc(sizeof(char)*256);
	char *subSampleName = malloc(sizeof(char)*256);
	int nSNP, nIND, K = 2, runSVD = 0, tmp, nF, nfile;
	int na, err = 0, sc = 1, i;
	double *Genotypes, *U, *V, *Sigma, *SNPSd, *Cov, *mAF, *miss, prop = .01; 
	int nSNP_file[MAXFILE];

	/* handleparams: gestion des paramètres */
	if (argc > 2){Welcome__f(0); if(handleParams__f(argc, argv, &nSNP, &nIND, &Genotypes, GenoFileName, &OutputFileName, &K, &runSVD, &subSampleName, &sc, &nfile, nSNP_file, &prop)) return 0;} else {Welcome__f(1); return 0;}

	nF = nIND > nSNP ? nSNP : nIND;
	nF = K;
	/* Allocation de mémoire */
        initializeVariables__f(&U, &Sigma, &V, &SNPSd, &Cov, &miss, &mAF, nF, nSNP, nIND);
        if(!strcmp(OutputFileName, "")) OutputFileName = "FastPCAdapt_output";
        printf("results in files %s*\n", OutputFileName);

	/* Calcule la matrice de covariance nxn */
	err = Cov_line(Cov, SNPSd, nSNP, nSNP_file, nIND, sc, GenoFileName, nfile);
	if (err) return 0;
printf("Covariance estimated\n");

	/* diagonalisation de la matrice nxn */
	/* renvoie la racine des valeurs propres dans Sigma */
	/* renvoie la matrice V (n x K) */
	diagonalize(Cov, nIND, K, Sigma, V);
/* V is a nxK matrix */
	/* calcule la matrice U (K x p) */
	err = Load_line(U, Sigma, V, miss, mAF, nSNP, nSNP_file, K, nIND, sc, GenoFileName, nfile);
	if (err) return 0;
printf("loadings computed\n");
	/* now we refer to V as tV, the K * nIND matrix */
	tr(V, nIND, K);

	/* Calcule les stats de l'ACP sur les loadings */
        getstat(U, Sigma, V, SNPSd, miss, mAF, nSNP, K, nIND, nF, sc, nSNP_file, GenoFileName, OutputFileName, nfile, prop);
        writeMatrix__f(U, Sigma, V, nSNP, K, nIND, nF, OutputFileName);

	free(miss);
	free(mAF);
	free(Cov);
	return 0;

}
