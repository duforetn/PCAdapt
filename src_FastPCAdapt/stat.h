
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

/* compare
 * generic function used in mad.
 */
int compare (const void* p1, const void* p2);

/* mad
 * Calculate the median absolute deviation of X, and sort X if reorder.
 */
double mad(double *X, int n, int reorder);

/* colSdev
 * Calculate standarddeviations of the colummns of Matrix in colSd.
 */
void colSdev(double *Matrix, int nrow, int ncol, double *colSd);

/* colMad
 * Calculate median absolute deviations of the columns of Matrix in colMad.
 */
void colMad(double *Matrix, int nrow, int ncol, double *colMad);

/* statmax
 * Calculate for each SNP the max of statndardized squarred loadings (in U). 
 */
void statmax(double *U, double *Sigma, double *Stats, double *colSd, int nSNP, int K, int nF, int nstat, int stat);

/* statsummad
 * Calculate for each SNP the sum of standardized squarred loadin
gs (in U). h'
 */
void statsummad(double *U, double *Stats, double *colSd, int nSNP, int K, int nF, int nstat, int stat);

/* statR
 * Calculate for each SNP the sum of standardized squarred loadin
gs (in U) times the eigenvalues (in Sigma). h
 */
void statR(double *U, double *Sigma, double *SNPSd, double *Stats, int nSNP, int nIND, int K, int nF, int nstat, int stat, int sc);

/* getstat
 * Calculate all the statistics for each SNPs, and write the results in files fileName(.*).
 */
void getstat(double *U, double *Sigma, double *V, double *SNPSd, double *miss, double *mAF, int nSNP, int K, int nIND, int nF, int sc, int nSNP_file[MAXFILE], char **GenoFileName, char *fileName, int nfile, double prop);

