#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../src_Lapack/lpk.h"

/* invMatrix
 * Calculate the inverse of Matrix, in Matrix.
 */
void invMatrix(double *Matrix, int nrow);

/* invSPDMatrix
 * Calculate the inverse of Matrix, in Matrix, if Matrix is a symmetric positive definite matrix.
 */
void invSPDMatrix(double *Matrix, int nrow);

/* logmatrixDet
 * Calculate the log of the determinant of Matrix.
 */
double logmatrixDet(double *Matrix, int nrow);

/* cholesky
 * Calculate Cholesky decomposition of A, with Lower matrix result in L.
 */
void cholesky(double* A, int n, double* L, int normalize);

/* SVD
 * Calculate Singular value decomposition of G, and writes U in Lambda, and Sigma * V in Factors.
 */
double SVD(double *G, double *Factors, double *Lambda, int K, int nSNP, int nIND);

/* tAA
 * Optimized way to calculate A tranpose * A in tAA.
 */
void tAA(double *A, double *tAA, int nrow, int ncol);

/* covariance
 * Calculate the covariance of A in Sigma, using tAA. If transpose == 0, calculate A * A transpose.
 */
void covariance(double *A, double *Sigma, int nrow, int ncol, int transpose);

/* diagonalize
 * Diagoanlize matrix cov and write squarred root eigenvalues in val, and eigenvectors in vect.
 */
void diagonalize(double *cov, int N, int K, double *val, double *vect);

