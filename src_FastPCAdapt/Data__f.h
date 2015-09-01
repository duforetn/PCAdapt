#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <unistd.h>
#include <string.h>

#define PRECISION_FLOAT 16
#define NA 9
#define MAXFILE 64

void Welcome__f(int help);

int initializeVariables__f(double **U, double **Sigma, double **V, double **SNPSd, double **Cov, double **miss, double **mAF, int K, int nSNP, int nIND);

int nboffile(char *argv[], int argc, int i, char *GenoFileName[]);

int handleParams__f(int argc, char *argv[], int *nSNP, int *nIND, double **Genotypes, char **GenoFileName, char **OutputFileName, int *K, int *runSVD, char **subSampleName, int *sc, int *nf, int *nSNP_file, double *prop);

void writeStats_multifile(char **GenoFileName, char *fileName, double *matrix, double *mAF, double *SNPSd, double *miss, int nrow, int nstat, int nfile, int nSNP_file[MAXFILE], double threshold);

void writeStats(char *fileName, double *matrix, double *mAF, double *SNPSd, double *miss, int nrow, int nstat, double threshold);

void writeMatrix__f(double *U, double *Sigma, double *V, int nSNP, int K, int nIND, int nF, char *OutputFileName);

