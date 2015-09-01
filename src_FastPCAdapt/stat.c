/*
   FastPCAdapt stat.c
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "stat.h"
#include "../src_Lapack/lpk.h"
#include "../src/linAlgebra.h"
#include "../src/matrix.h"
#include "../src/Data.h"
#include "Data__f.h"

 	
int compare (const void* p1, const void* p2){

	if ( *((double *) p1) < *((double *) p2)) return -1;
	if ( *((double *) p1) == *((double *) p2)) return 0;
	if ( *((double *) p1) > *((double *) p2)) return 1;

}

double mad(double *X, int n, int reorder){

	double medX, mad, tmp;
	int i, imin, start;
	if (reorder){
		qsort((void *) X, (size_t) n, sizeof(double), compare);
		if (n%2) medX = X[(n - 1)/2];
		else medX = X[n/2]/2.0 + X[n/2 - 1]/2.0;
		for (i=0; i<n; i++) X[i] = fabs(X[i] - medX);
		qsort((void *) X, (size_t) n, sizeof(double), compare);
                if (n%2) mad = X[n/2];
                else mad = (X[n/2] + X[n/2 + 1])/2.0;
	} else {
		double *myX = malloc(sizeof(double)*n);
		for (start=0; start<n - 1; start++){
			imin = start;
			for (i=start; i<n; i++){
				if (X[imin] > X[i]) imin = i;
			}
			myX[start] = X[imin];
		}
                if (n%2) medX = myX[n/2];
                else medX = (myX[n/2] + myX[n/2 + 1])/2;
                for (i=0; i<n; i++) myX[i] = fabs(myX[i] - medX);
                qsort((void *) myX, (size_t) n, sizeof(double), compare);
                if (n%2) mad = myX[n/2];
                else mad = (myX[n/2] + myX[n/2 + 1])/2;
		free(myX);
	}
	return mad;
}

void colSdev(double *Matrix, int nrow, int ncol, double *colSd){

        double *cMeans = malloc(sizeof(double)*ncol);
        int i, j;
        colMeans(Matrix, nrow, ncol, cMeans);

        for (i=0; i<ncol; i++){
		colSd[i] = 0;
                for (j=0; j<nrow; j++){
                        colSd[i] += (Matrix[i*ncol + j] - cMeans[i])*(Matrix[i*ncol + j] - cMeans[i]);
                }
                colSd[i] = sqrt(colSd[i]/(nrow - 1));
        }
        free(cMeans);
}

void colMad(double *Matrix, int nrow, int ncol, double *colMad){

        int i, j;
        double *col = malloc(sizeof(double)*nrow);
        for (i=0; i<ncol; i++){
                for (j=0; j<nrow; j++){
                        col[j] = Matrix[j*ncol + i];
                }
                colMad[i] = mad(col, nrow, 1);
        }
        free(col);
}

void statmax(double *U, double *Sigma, double *Stats, double *colSd, int nSNP, int K, int nF, int nstat, int stat){

        int i, j, jmax;
        double S, Smax;
        double scale = (double) sqrt((double) nSNP);
        for (i=0; i<nSNP; i++){
                jmax = 0;
                Smax = 0;
                for (j=0; j<K; j++){
                        S = U[i*nF + j]*U[i*nF + j]*scale*scale/(colSd[j]*colSd[j]);
                        if (S > Smax){
                                jmax = j;
                                Smax = S;
                        }
                }
                Stats[i*nstat + stat] = Smax;
                Stats[i*nstat + stat + 1] = jmax + 1;
        }

}

void statsummad(double *U, double *Stats, double *colSd, int nSNP, int K, int nF, int nstat, int stat){

        int i, j;
	double S;
	double scale = (double) sqrt((double) nSNP);
        for (i=0; i<nSNP; i++){
		S = 0;
		for (j=0; j<K; j++){
			S += (U[i*nF + j]*scale)*(U[i*nF + j]*scale)/(colSd[j]*colSd[j]);
		}
		Stats[i*nstat + stat] = S;
        }
}

void statR(double *U, double *Sigma, double *SNPSd, double *Stats, int nSNP, int nIND, int K, int nF, int nstat, int stat, int sc){

        int k, j;
        double d;
        for (j=0; j<nSNP; j++){
                Stats[j*nstat + stat] = 0;
                d = 0;
                for (k=0; k<K; k++){
                        Stats[j*nstat + stat] += Sigma[k]*Sigma[k]*U[j*nF + k]*U[j*nF + k]/(nIND - 1);
                }
                if (!sc){
			if (SNPSd[j] == 0){
				 Stats[j*nstat + stat] = 0;
			} else {
				Stats[j*nstat + stat] = sqrt(Stats[j*nstat + stat]/(SNPSd[j]*SNPSd[j]));
			}
		}
        }

}

void getstat(double *U, double *Sigma, double *V, double *SNPSd, double *miss, double *mAF, int nSNP, int K, int nIND, int nF, int sc, int nSNP_file[MAXFILE], char **GenoFileName, char *fileName, int nfile, double prop){

        int nstat = 4, i;
	double threshold;
        double *Stats = malloc(sizeof(double)*nSNP*(nstat + 2));
	double *Stat_cpy = malloc(sizeof(double)*nSNP);
        double *colSd = calloc(nF, sizeof(double));
//        colMad(U, nSNP, nF, colSd);
//        colSdev(U, nSNP, nF, colSd);
//displayMatrix(colSd, 1, nF);
        for (i=0; i<nF; i++) if(colSd[i] == 0) colSd[i] = 1.0;

        statsummad(U, Stats, colSd, nSNP, K, nF, nstat, 0);
        statmax(U, Sigma, Stats, colSd, nSNP, K, nF, nstat, 1);
        statR(U, Sigma, SNPSd, Stats, nSNP, nIND, K, nF, nstat, 3, sc);

	for (i=0; i<nSNP; i++) Stat_cpy[i] = Stats[i*nstat];	
        qsort((void *) Stat_cpy, (size_t) nSNP, sizeof(double), compare);
        threshold = Stat_cpy[(int) floor(nSNP*(1 - prop))];
        if (nfile == 1) {writeStats(fileName, Stats, mAF, SNPSd, miss, nSNP, nstat, threshold);}
	else {writeStats_multifile(GenoFileName, fileName, Stats, mAF, SNPSd, miss, nSNP, nstat, nfile, nSNP_file, threshold);}
	free(Stat_cpy);
        free(colSd);
        free(Stats);
}


