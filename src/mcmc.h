#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "matrix.h"
#include "random.h"

/* initializeVariables
 * Initialize and allocate memory for the matrices and vectors of the MCMC.
 */
int initializeVariables(double **Factors, double **sigma2_F, double **Lambda, double **Pz, double **Piz, double **C2, double **sigma2_Lambda, double **rowSd, int K, double alpha, int nSNP, int nIND);

/* smatrInit
 * Read data in file subsampleName, and perform SVD. U is store in Lambda, and eigenvalues in sigma2_F.
 */
double smartInit(double *Lambda, double *sigma2_F, int nIND, int K, char *subSampleName);

/* runMCMC
 * run the MCMC algorithm on allocated data, and writes results.
 */
int runMCMC(double *G, int *missing, double *Factors, double *sigma2_F, double *Lambda, double *Pz, double *Piz, double *sigma2, int K, double alpha, double Tau2, double *C2, double *sigma2_Lambda, double Isingb, int nSNP, int nIND, int nsteps, int burnin, int na, int prior, int runSVD, char *OutputFileName, char *subSampleName, double prop);

