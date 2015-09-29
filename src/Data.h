#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <unistd.h>
#include <string.h>

#define PRECISION_FLOAT 16
#define NA 9

/* Welcome
 *Display help menu
 */
void Welcome(int help);

/* extractData
 * Store data from file NameINPUT in Matrix. format used is "genotype".
 */
int extractData(char *NameINPUT, int nrow, int ncol, double *Matrix, char *format);

/* DataImputation
 * Replace missing data in G. method "naive" draws allele from allele frequencies
 */
int DataImputation(double *G, int **missing, int nrow, int ncol, char *method);

/* DataCounter
 * Count the number of variables and individuals in NameINPUT.
 */
int DataCounter(char *NameINPUT, int *nrow, int *ncol, char *format);

/* handleParams
 * Deal with parameters in the command line as described in Help.
 */
int handleParams(int argc, char *argv[], int *nSNP, int *nIND, double **Genotypes, char **GenoFileName, char **OutputFileName, int *K, int *nsteps, int *burnin, int *runSVD, char **subSampleName, int *transpose, int *prior, double *b, int *sc, double *prop);

/* getqval
 * Estimate q-values from a nSNP vector of p-values Proba.
 */
void getqval(double *Proba, double *qval, int nSNP, int K);

/* writeResults
 * Write scores for each snps in fileName, and parameter estimates from the MCMC.
 */
void writeResults(char *fileName, double *matrix, double *BF, double *PO, int nSNP, int K, double err, double alpha, double* C2, double *sigma2_F, int *order_v);

/* writeScores
 * Write U matrix in fileName.scores.
 */
void writeScores(char *fileName, double *matrix, double *loadings, int K, int nIND, int nSNP, int *order_v);

/* compareBF
 * Generic function used in getqval.
 */
int compareBF(const void* p1, const void* p2);

/* order
 * Sort vector sigma2_F in decreaasing order.
 */
void order(double *sigma2_F, int *order_v, int K);

/* writeBF
 * Get the top prop Bayes factors, and write the result in fileName.topBF.
 */
void writeBF(char *fileName, double *BF, double *Pz_chain, int nSNP, int K, double prop, int *order_v);
