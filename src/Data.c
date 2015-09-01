/*
   PCAdapt Data.c
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
#include "Data.h"
#include "random.h"

#define PRECISION_FLOAT 16
#define NA 9

void Welcome(int help){
	printf("\t\t/************************\\\n");
	printf("\t\t|*** PCAdapt software ***|\n");
	printf("\t\t\\************************/\n\n");
	printf("\n\n");

	if (help) printf("Command line Help:\n\t-i path and name of the file containing your genotypes.\n\t-K number of factors to use\n\t-b number of burnin iterations\n\t-s number of iteration (including burnin)\n\t-o name of the output file\n\t-t 1 if you have individuals in lines\n\t-B Ising model parameter\n\t-I name of data file to initialise the chain with SVD\n\t-S 0 if you want the data to be unscaled.\n");
}

int extractData(char *NameINPUT, int nrow, int ncol, double *Matrix, char *format){

        int i = 0, j = 0;
	float value;
	FILE *INPUTFile;
        if((INPUTFile = fopen(NameINPUT, "r")) == NULL){
                printf("Error, invalid input file\n");
                return 1;
        }

	if(!strcmp(format, "genotypes")){
                while((fscanf(INPUTFile, "%g", &value) != EOF)){
                        Matrix[i*ncol + j] = (double) value;
	//		printf("%i, %i -> %g\n", i, j, value);
                        j++;
                        if(j == ncol){
                                j = 0;
                                i++;
                        }
                }
	}
	fclose(INPUTFile);
	return 0;
}

int DataImputation(double *G, int **missing, int nrow, int ncol, char *method){

        int i, j, geno, na = 0, count = 0;
        double r;
        double *AF = calloc(nrow, sizeof(double));
	for (i=0; i<ncol*nrow; i++) if (G[i] == NA) na++;
        *missing = calloc(2*na, sizeof(int));

        if (!strcmp(method, "naive")){
                rowMeans(G, nrow, ncol, AF);
                for (i=0; i<nrow; i++){
                        for (j=0; j<ncol; j++){
                                geno = 0;
                                if (G[i*ncol + j] == NA){
                                        r = drand();
                                        if(r < AF[i]/2) geno++;
                                        r = rand_double(0, 1);
                                        if(r < AF[i]/2) geno++;
                                        G[i*ncol + j] = geno;
	        	                (*missing)[count*2] = i;
		                        (*missing)[count*2 + 1] = j;
					count++;
                                }
                        }
                }
        }
        free(AF);
        return na;
}


int DataCounter(char *NameINPUT, int *nrow, int *ncol, char *format){

        FILE *INPUTFile;
        long long  nbbn = 0;
        long long nbsp = 0;
        long long res;
        int currentchar;
        int prevchar;
        if((INPUTFile = fopen(NameINPUT, "r")) == NULL){
                printf("Error, invalid input file\n");
                return 1;
        }
        currentchar = fgetc(INPUTFile);
        while(currentchar != EOF){
                if (currentchar == 10){
                        nbbn++;
                        if (prevchar != 32 && prevchar != '\t'){
                                nbsp++;
                        }
                }
                if ((currentchar == 32 || prevchar == '\t') && (prevchar != 32 || prevchar != '\t')) nbsp++;
                prevchar = currentchar;
                currentchar = fgetc(INPUTFile);
        }

        fclose(INPUTFile);
        if(!strcmp(format, "genotypes")){
                *nrow = (int)  nbbn;
                res = nbsp/nbbn;
                *ncol = res;
        }
        return 0;

}

int handleParams(int argc, char *argv[], int *nSNP, int *nIND, double **Genotypes, char **GenoFileName, char **OutputFileName, int *K, int *nsteps, int *burnin, int *runSVD, char **subSampleName, int *transpose, int *prior, double *b, int *sc, double *prop){

	int i, tmp;
        for (i=0; i<argc; i++){
                if (argv[i][0] == '-'){
                        switch(argv[i][1]){
	                        case 'i':
        	                        *GenoFileName = strdup(argv[i + 1]);
                	                printf("Genotypes in file %s\n", *GenoFileName);
                        	        if(!DataCounter(*GenoFileName, nSNP, nIND, "genotypes")){
//                              	          printf("%i snps typed on %i individuals\n", *nSNP, *nIND);
                                        	*Genotypes = calloc((*nSNP)*(*nIND), sizeof(double));

						int res;
        	                                res = extractData(*GenoFileName, *nSNP, *nIND, *Genotypes, "genotypes");
                	                } else {return 1;}
	                        break;
				case 'K':
					*K = atoi(argv[i + 1]);
				break;
				case 's':
					*nsteps = atoi(argv[i + 1]);
				break;
				case 'b':
					*burnin = atoi(argv[i + 1]);
				break;
				case 'o':
					*OutputFileName = strdup(argv[i + 1]);
				break;
				case 'I':
					*runSVD = 1;
					*subSampleName = strdup(argv[i + 1]);
				break;
				case 't':
					*transpose = atoi(argv[i + 1]);
				break;
        	                case 'p':
                        	        *prop = (double) atof(argv[i + 1]);
	                        break;
				case 'B':
					*b = (double) atof(argv[i + 1]);
				break;
				case 'S':
					*sc = atoi(argv[i + 1]);
				break;
			}
		}
        }
        if (*transpose == 1){
                tr(*Genotypes, *nSNP, *nIND);
                tmp = *nSNP;
                *nSNP = *nIND;
                *nIND = tmp;
                printf("%i snps typed on %i individuals\n", *nSNP, *nIND);
        } else {
                printf("%i snps typed on %i individuals\n", *nSNP, *nIND);
        }

	if (*K >= *nSNP || *K >= *nIND){
		printf("***WARNING*** K must be lower than the rank of the design matrix\n");
		*K = *nSNP < *nIND ? *nSNP - 1 : *nIND - 1;
	}
	return 0;
}

void getqval(double *Proba, double *qval, int nSNP, int K){

	int i, snp, imin = 0;
	double q = 0;
	double *Pz = calloc(nSNP, sizeof(double));
	for (i=0; i<nSNP; i++) Pz[i] = Proba[i*(K + 1)];
	for (snp=0; snp<nSNP; snp++){
		for (i=0; i<nSNP; i++){
			if (Pz[i] < Pz[imin]) imin = i;
		}	
		q = (q*snp + Pz[imin])/(snp + 1);
		qval[imin] = q;
		Pz[imin] = 12;
	}
	free(Pz);
}

void writeResults(char *fileName, double *matrix, double *BF, double *PO, int nSNP, int K, double err, double alpha, double* C2, double *sigma2_F, int *order_v){

	int i, j;
	int nrow = nSNP, ncol = K + 1;
	double *qval = calloc(nSNP, sizeof(double));
	FILE *resultFile;
        char *fileName_ext = calloc(256, sizeof(char));
        strcpy(fileName_ext, fileName);
        strcat(fileName_ext, ".stats");
	printf("Opening file %s...\n", fileName_ext);
	if((resultFile = fopen(fileName_ext, "w")) == NULL) printf("ERROR, unable to open %s\n", fileName_ext);
        fprintf(resultFile, "error: %f\n", err);
	fprintf(resultFile, "pi: %f\n", alpha);
	fprintf(resultFile, "c2: ");
	for (i=0; i<K; i++) fprintf(resultFile, "%f ", C2[order_v[i]]);
	fprintf(resultFile, "\n");
	
	fprintf(resultFile, "rho2:\n");

	/* get q values */
	//getqval(matrix, qval, nSNP, K);

        for (i=0; i<ncol - 1; i++) fprintf(resultFile, "%lf ", sigma2_F[order_v[i]]);

	fclose(resultFile);

        printf("Opening file %s...\n", fileName);
        if((resultFile = fopen(fileName, "w")) == NULL) printf("ERROR, unable to open %s\n", fileName);
	fprintf(resultFile, "logBF\tlogPO\t");
	for (i=1; i<ncol; i++) fprintf(resultFile, "P(Z = %d|O)\t", i);
//	fprintf(resultFile, "q-value\n");
        fprintf(resultFile, "\n");
	for (i=0; i<nrow; i++){
		double c = .0;
		for (j=1; j<ncol; j++) c += matrix[i*ncol + j];
		fprintf(resultFile, "%f\t%f\t", BF[i], PO[i]);
		for (j=1; j<ncol; j++){
			fprintf(resultFile, "%f\t", matrix[i*ncol + order_v[j - 1] + 1]/c);
		}
//		fprintf(resultFile, "%f\n", qval[i]);
                fprintf(resultFile, "\n");
	}
	fclose(resultFile);
	free(qval);
}

void writeScores(char *fileName, double *matrix, double *loadings, int K, int nIND, int nSNP, int *order_v){
	FILE *resultFile, *loadingsFile;
	int i, j;
	char *fileName_ext = calloc(256, sizeof(char));
        char *fileName_loadings = calloc(256, sizeof(char));
	strcpy(fileName_ext, fileName);
	strcat(fileName_ext, ".scores");
        printf("Opening file %s...\n", fileName_ext);
        if((resultFile = fopen(fileName_ext, "w")) == NULL) printf("ERROR, unable to open %s\n", fileName_ext);
	for (i=0; i<K; i++){
		for (j=0; j<nIND; j++){
			fprintf(resultFile, "%f ", matrix[order_v[i]*nIND + j]);
		}
		fprintf(resultFile, "\n");
	}
        strcpy(fileName_loadings, fileName);
        strcat(fileName_loadings, ".loadings");
        printf("Opening file %s...\n", fileName_loadings);
        if((loadingsFile = fopen(fileName_loadings, "w")) == NULL) printf("ERROR, unable to open %s\n", fileName_loadings);
        for (i=0; i<nSNP; i++){
                for (j=0; j<K; j++){
                        fprintf(loadingsFile, "%f ", loadings[i*K + order_v[j]]);
                }
                fprintf(loadingsFile, "\n");
        }

	
}

int compareBF(const void* p1, const void* p2){

        if ( *((double *) p1) < *((double *) p2)) return -1;
        if ( *((double *) p1) == *((double *) p2)) return 0;
        if ( *((double *) p1) > *((double *) p2)) return 1;

}

void order(double *sigma2_F, int *order_v, int K){

	int i, imax, start = 0;
	double max;
	int *tab = calloc(K, sizeof(int));
	for (start=0; start<K; start++){
		max = -1;
		for (i=0; i<K; i++){
			if( tab[i] == 0 && sigma2_F[i] > max){
				max = sigma2_F[i];
				imax = i;
			}
		}
		tab[imax] = 1;
		order_v[start] = imax;
	}

}

void writeBF(char *fileName, double *BF, double *Pz_chain, int nSNP, int K, double prop, int *order_v){

        int i, j, Factor, rank;
	double threshold;
        FILE *resultFileBF;
        char *fileName_BF = calloc(256, sizeof(char));
        strcpy(fileName_BF, fileName);
        strcat(fileName_BF, ".topBF");
	double *BF_cpy = malloc(sizeof(double)*nSNP);
	for (i=0; i<nSNP; i++) BF_cpy[i] = BF[i];

	qsort((void *) BF_cpy, (size_t) nSNP, sizeof(double), compareBF);
	threshold = BF_cpy[(int) floor(nSNP*(1 - prop))];

        printf("Opening file %s...\n", fileName_BF);
        if((resultFileBF = fopen(fileName_BF, "w")) == NULL) printf("ERROR, unable to open %s\n", fileName_BF);

	fprintf(resultFileBF, "snp\tFactor\tBF\n");

	for (i=0; i<nSNP; i++){
		Factor = findMax(Pz_chain + i*(K + 1) + 1, 1, K);
		j = 0;
		while(Factor != order_v[j]) j++;
		if (BF[i] >= threshold){
			fprintf(resultFileBF, "%i\t%i\t%f\n", i + 1, j + 1, pow(10, BF[i]));
		}
	}
	free(BF_cpy);
        fclose(resultFileBF);
}

