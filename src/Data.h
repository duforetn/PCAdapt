#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <unistd.h>
#include <string.h>

#define PRECISION_FLOAT 16
#define NA 9

void Welcome(int help);

int extractData(char *NameINPUT, int nrow, int ncol, double *Matrix, char *format);

int DataImputation(double *G, int **missing, int nrow, int ncol, char *method);

int DataCounter(char *NameINPUT, int *nrow, int *ncol, char *format);

int handleParams(int argc, char *argv[], int *nSNP, int *nIND, double **Genotypes, char **GenoFileName, char **OutputFileName, int *K, int *nsteps, int *burnin, int *runSVD, char **subSampleName, int *transpose, int *prior, double *b, int *sc, double *prop);

void getqval(double *Proba, double *qval, int nSNP, int K);

void writeResults(char *fileName, double *matrix, double *BF, double *PO, int nSNP, int K, double err, double alpha, double* C2, double *sigma2_F, int *order_v);

void writeScores(char *fileName, double *matrix, double *loadings, int K, int nIND, int nSNP, int *order_v);

int compareBF(const void* p1, const void* p2);

void order(double *sigma2_F, int *order_v, int K);

void writeBF(char *fileName, double *BF, double *Pz_chain, int nSNP, int K, double prop, int *order_v);
