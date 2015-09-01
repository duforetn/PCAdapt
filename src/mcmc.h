#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "matrix.h"
#include "random.h"

int initializeVariables(double **Factors, double **sigma2_F, double **Lambda, double **Pz, double **Piz, double **C2, double **sigma2_Lambda, double **rowSd, int K, double alpha, int nSNP, int nIND);

double smartInit(double *Lambda, double *sigma2_F, int nIND, int K, char *subSampleName);

int runMCMC(double *G, int *missing, double *Factors, double *sigma2_F, double *Lambda, double *Pz, double *Piz, double *sigma2, int K, double alpha, double Tau2, double *C2, double *sigma2_Lambda, double Isingb, int nSNP, int nIND, int nsteps, int burnin, int na, int prior, int runSVD, char *OutputFileName, char *subSampleName, double prop);

