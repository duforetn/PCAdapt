#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../src_Lapack/lpk.h"

void invMatrix(double *Matrix, int nrow);

void invSPDMatrix(double *Matrix, int nrow);

double logmatrixDet(double *Matrix, int nrow);

void cholesky(double* A, int n, double* L, int normalize);

double SVD(double *G, double *Factors, double *Lambda, int K, int nSNP, int nIND);

void tAA(double *A, double *tAA, int nrow, int ncol);

void covariance(double *A, double *Sigma, int nrow, int ncol, int transpose);

void diagonalize(double *cov, int N, int K, double *val, double *vect);

