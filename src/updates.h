#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "matrix.h"
#include "random.h"


void updateMu(double *mu, double *G, double *Factors, double *Lambda, double sigma2, int nSNP, int nIND, int K);

void updateLambda(double *Lambda, double *G, double *Factors, double sigma2, double *sigma2_Lambda, int K, int nSNP, int nIND);

void updatePiz(double *Piz, int *Z, double alpha, int K, int nSNP);

double ConstantPriorOdds();

void updatePz(double *Pz, double *Factors, double alpha, double *sigma2_F, int *Z, double Tau2, double *C2, double sigma2, double Isingb, int nSNP, int K, int prior, double *PriorOdds, double *PO, int step, int burnin);

int updateZ(int *Z, double *Pz, int nSNP, int K);

void updateAlpha(double *alpha, int noutlier, int nSNP);

void updateFactors(double *Factors, int *Z, double *Lambda, double *G, double Tau2, double *C2, double *sigma2_F, double sigma2, int nSNP, int nIND, int K, int prior);

int updateSigma2_F(double *sigma2_F, double *Factors, double sigma2, int *Z, double Tau2, double *C2, int K, int nSNP, int prior);

void updateSigma2_Lambda(double *sigma2_Lambda, double *Lambda, int K, int nIND);

void updateSigma2(double *sigma2, double *G, double *Factors, double *Lambda, double *sigma2_F, double *C2, int *Z, int nSNP, int nIND, int K, int prior);

void updateTau2(double *Tau2, double *C2, double *sigma2_F, double *Factors, double sigma2, int *Z, int K, int nSNP, int prior);
