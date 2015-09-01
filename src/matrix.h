#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../src_Lapack/lpk.h"

void displayMatrix(double *Matrix, int nrow, int ncol);

void rowMeans(double *Matrix, int nrow, int ncol, double *Means);

void colMeans(double *Matrix, int nrow, int ncol, double *Means);

void sumMatrix(double *A, double *B, double *res, int nrow, int ncol);

void difMatrix(double *A, double *B, double *res, int nrow, int ncol);

void prodMatrix(double *A, double *B, double *res, int nrowA, int ncolA, int nrowB, int ncolB);

void tr(double *Matrix, int nrow, int ncol);

int findMax(double *M, int nrow, int ncol);

int findMin(double *M, int nrow, int ncol);

/* Centers and scales (if std) a matrix by lines */
void scale(double *Matrix, double *rowSd, int nrow, int ncol, int center, int std);

