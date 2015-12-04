
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../src_Lapack/lpk.h"
#include "../src/linAlgebra.h"
#include "../src/matrix.h"
#include "../src/Data.h"
#include "Data__f.h"
#define NA 9
#define MAXFILE 64

int compare (const void* p1, const void* p2);

double mad(double *X, int n, int reorder);

void colSdev(double *Matrix, int nrow, int ncol, double *colSd);

void colMad(double *Matrix, int nrow, int ncol, double *colMad);

void statmax(double *U, double *Simga, double *Stats, double *colSd, int nSNP, int K, int nF, int nstat, int stat);

void statsummad(double *U, double *Stats, double *colSd, int nSNP, int K, int nF, int nstat, int stat);

void statcor(double *Genotypes, double *Stats, double *Sigma, double *V, int nSNP, int K, int nIND, int nF, int nstat, int stat);

void statR(double *U, double *Sigma, double *SNPSd, double *Stats, int nSNP, int nIND, int K, int nF, int nstat, int stat, int sc);

void getstat(double *U, double *Sigma, double *V, double *SNPSd, double *miss, double *mAF, int nSNP, int K, int nIND, int nF, int sc, int nSNP_file[MAXFILE], char **GenoFileName, char *fileName, int nfile, double prop);

