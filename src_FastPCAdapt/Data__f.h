#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <unistd.h>
#include <string.h>

#define PRECISION_FLOAT 16
#define NA 9
#define MAXFILE 64

/* Welcome
 *Display help menu of the fast version.
 */
void Welcome__f(int help);

/* initializeVariables__f
 * Initialize and allocate memory for the matrices and vectors of the fast version.
 */
int initializeVariables__f(double **U, double **Sigma, double **V, double **SNPSd, double **Cov, double **miss, double **mAF, int K, int nSNP, int nIND);

/* nboffile
 * Calculate the number of input files in the command line, and return it.
 */
int nboffile(char *argv[], int argc, int i, char *GenoFileName[]);

/* handleParams__f
 * Deal with parameters in the command line as described in Help, for the fast version.
 */
int handleParams__f(int argc, char *argv[], int *nSNP, int *nIND, double **Genotypes, char **GenoFileName, char **OutputFileName, int *K, int *runSVD, char **subSampleName, int *sc, int *nf, int *nSNP_file, double *prop);

/* writeStats_multifile
 * Multifile version of writeStats. one file per inputfile is written.
 * Write the nstat statistics for each snp in fileName(.topBF). threshold is the specified threshold for the .topSNP file. 
 */

void writeStats_multifile(char **GenoFileName, char *fileName, double *matrix, double *mAF, double *SNPSd, double *miss, int nrow, int nstat, int nfile, int nSNP_file[MAXFILE], double threshold);

/* writeStats
 * Write the nstat statistics for each snp in fileName(.topBF). threshold is the specified threshold for the .topSNP file. 
 */
void writeStats(char *fileName, double *matrix, double *mAF, double *SNPSd, double *miss, int nrow, int nstat, double threshold);

/* writeMatrix__f
 * Write matrices U, SIgma and V in files OutputFileName.(scores, sigma, loadings).
 */
void writeMatrix__f(double *U, double *Sigma, double *V, int nSNP, int K, int nIND, int nF, char *OutputFileName);

