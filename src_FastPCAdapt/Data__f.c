/*
   FastPCAdapt Data__f.c
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
#include "Data__f.h"
#include "../src/random.h"

#define PRECISION_FLOAT 16
#define NA 9

void Welcome__f(int help){
	printf("\t\t/****************************\\\n");
	printf("\t\t|*** FastPCAdapt software ***|\n");
	printf("\t\t\\****************************/\n\n");
	printf("\n\n");

	if (help) printf("Command line Help:\n\t-i paths and names of the files containing your genotypes, it can be several files.\n\t-K number of factors to use\n\t-o name of the output file\n\t-S 0 do not scale the snps\n\n");
}

int initializeVariables__f(double **U, double **Sigma, double **V, double **SNPSd, double **Cov, double **miss, double **mAF, int K, int nSNP, int nIND){

        int i, j;
        init_random();

	*Sigma = calloc(K, sizeof(double));
        *U = malloc(sizeof(double)*nSNP*K);
        for (i=0; i<K; i++){
                for (j=0; j<nSNP; j++){
                        *(*(U) + i*nSNP + j) = rand_normal(0, 1);
                }
        }
	*SNPSd = malloc(sizeof(double)*nSNP);
        *V = calloc(K*nIND, sizeof(double));
	*Cov = calloc(nIND*nIND, sizeof(double));
        for (i=0; i<K*nIND; i++) *(*(V) + i) = rand_normal_r();
	*miss = calloc(nSNP, sizeof(double));
	*mAF = calloc(nSNP, sizeof(double));

        return 0;
}

int nboffile(char *argv[], int argc, int i, char *GenoFileName[]){

	int n = 0;
	while ((n + i < argc) && (argv[n + i][0] != '-')){
		GenoFileName[n] = strdup(argv[n + i]);
		n++;
	}
	return n;
}

int handleParams__f(int argc, char *argv[], int *nSNP, int *nIND, double **Genotypes, char **GenoFileName, char **OutputFileName, int *K, int *runSVD, char **subSampleName, int *sc, int *nf, int *nSNP_file, double *prop){

	int i, j, tmp, nfile, snp;
	*nSNP = 0;
        for (i=0; i<argc; i++){
                if (argv[i][0] == '-'){
                        switch(argv[i][1]){
	                        case 'i':
					nfile = nboffile(argv, argc, i + 1, GenoFileName);
                	                for(j=0; j<nfile; j++){
						printf("Reading file %s... ", GenoFileName[j]);
                        		        if(!DataCounter(GenoFileName[j], &snp, nIND, "genotypes")){
                	              	        	printf("%i snps typed on %i individuals\n", snp, *nIND);
							nSNP_file[j] = snp;
							*nSNP += snp;
	                	                } else {return 1;}
					}
					*nf = nfile;
					if (nfile > 1) printf("A total of %d snps typed on %d individuals in the %d files\n", *nSNP, *nIND, nfile);
	                        break;
				case 'K':
					*K = atoi(argv[i + 1]);
				break;
				case 'o':
					*OutputFileName = strdup(argv[i + 1]);
				break;
				case 'I':
					*runSVD = 1;
					*subSampleName = strdup(argv[i + 1]);
				break;
				case 'S':
					*sc = atoi(argv[i + 1]);
				break;
				case'p':
					*prop = (double) atof(argv[i + 1]);
				break;
			}
		}
        }
	if (*K >= *nSNP || *K >= *nIND){
		printf("***WARNING*** K must be lower than the rank of the design matrix\n");
		*K = *nSNP < *nIND ? *nSNP - 1 : *nIND - 1;
	}
	return 0;
}

void writeStats_multifile(char **GenoFileName, char *fileName, double *matrix, double *mAF, double *SNPSd, double *miss, int nrow, int nstat, int nfile, int nSNP_file[MAXFILE], double threshold){

        FILE *resultFile;
        FILE *topFile;
        char *topfileName = calloc(256, sizeof(char));
        strcpy(topfileName, topfileName);
        strcat(topfileName, ".top");

        int i, j, file, start = 0;
	char *filestat = malloc(sizeof(char)*256);	
	for (file=0; file < nfile; file++){

	        strcpy(filestat, fileName);
		char *ext = malloc(sizeof(char)*10);
		sprintf(ext, "_file%d", file + 1);
 	        strcat(filestat, ext);
	        printf("Opening file %s...\n", filestat);
	        strcpy(topfileName, filestat);
	        strcat(topfileName, ".top");
                printf("Opening file %s...\n", topfileName);
	        if((resultFile = fopen(filestat, "w")) == NULL) printf("ERROR, unable to open %s\n", filestat);
		if((topFile = fopen(topfileName, "w")) == NULL) printf("ERROR, unable to open %s\n", topfileName);

	        fprintf(resultFile, "h\'\taxis\th\tmAF\tmiss\n");
		fprintf(topFile, "snp\tPC\tscore\n");
        	for (i=start; i<start + nSNP_file[file]; i++){
	                for (j=0; j<nstat; j++){
				//
        	                if (j != 1) fprintf(resultFile, "%f ", matrix[i*nstat + j]);
                	}
	                fprintf(resultFile, "%f %f\n", mAF[i], miss[i]);
	                if (matrix[i*nstat] > threshold) {
        	                fprintf(topFile, "%d %d %f\n", i - start + 1, (int) matrix[i*nstat + 2], matrix[i*nstat]);
	                }			
        	}
		start += nSNP_file[file];
		fclose(resultFile);
		fclose(topFile);
	}
}

void writeStats(char *fileName, double *matrix, double *mAF, double *SNPSd, double *miss, int nrow, int nstat, double threshold){

        FILE *resultFile;
	FILE *topFile;
        char *topfileName = calloc(256, sizeof(char));
        strcpy(topfileName, fileName);
        strcat(topfileName, ".top");

        int i, j;
        printf("Opening file %s...\n", fileName);
        printf("Opening file %s...\n", topfileName);
        if((resultFile = fopen(fileName, "w")) == NULL) printf("ERROR, unable to open %s\n", fileName);
	if((topFile = fopen(topfileName, "w")) == NULL) printf("ERROR, unable to open %s\n", topfileName);
        fprintf(resultFile, "h\'\taxis\th\tmAF\tmiss\n");
	fprintf(topFile, "snp\tPC\tscore\n");
        for (i=0; i<nrow; i++){
                for (j=0; j<nstat; j++){
			//Don't write the max, column #1
                        if (j != 1) fprintf(resultFile, "%f ", matrix[i*nstat + j]);
                }
                fprintf(resultFile, "%f %f\n", mAF[i], miss[i]);
		if (matrix[i*nstat] >= threshold) {
			fprintf(topFile, "%d %d %f\n", i + 1, (int) matrix[i*nstat + 2], matrix[i*nstat]);
		}
        }
	
}


void writeMatrix__f(double *U, double *Sigma, double *V, int nSNP, int K, int nIND, int nF, char *OutputFileName){

        FILE *resultFile;	
	char *fileU = malloc(sizeof(char)*256);
        char *fileS = malloc(sizeof(char)*256);
        char *fileV = malloc(sizeof(char)*256);
        int i, j;
        strcpy(fileU, OutputFileName);
        strcat(fileU, ".loadings");

        printf("Opening file %s...\n", fileU);
        if((resultFile = fopen(fileU, "w")) == NULL) printf("ERROR, unable to open %s\n", fileU);
        for (i=0; i<nSNP; i++){
                for (j=0; j<nF; j++){
                        fprintf(resultFile, "%f ", U[i*nF + j]);
                }
                fprintf(resultFile, "\n");
        }
	fclose(resultFile);

        strcpy(fileS, OutputFileName);
        strcat(fileS, ".sigma");

        printf("Opening file %s...\n", fileS);
        if((resultFile = fopen(fileS, "w")) == NULL) printf("ERROR, unable to open %s\n", fileS);
        for (i=0; i<nF; i++){
        	fprintf(resultFile, "%f ", Sigma[i]/sqrt(nIND - 1));
        }
	fprintf(resultFile, "\n");
        fclose(resultFile);

        strcpy(fileV, OutputFileName);
        strcat(fileV, ".scores");

        printf("Opening file %s...\n", fileV);
        if((resultFile = fopen(fileV, "w")) == NULL) printf("ERROR, unable to open %s\n", fileV);
        for (i=0; i<nF; i++){
                for (j=0; j<nIND; j++){
                        fprintf(resultFile, "%f ", V[i*nIND + j]);
                }
                fprintf(resultFile, "\n");
        }
        fclose(resultFile);
	free(fileU);
        free(fileS);
        free(fileV);

}

