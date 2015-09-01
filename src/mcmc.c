/*
   PCAdapt mcmc.c
   Copyright (C) 2014 Nicolas Duforet Frebourg

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define DEBUG 0
#include "random.h"
#include "Data.h"
#include "matrix.h"
#include "updates.h"

int initializeVariables(double **Factors, double **sigma2_F, double **Lambda, double **Pz, double **Piz, double **C2, double **sigma2_Lambda, double **rowSd, int K, double alpha, int nSNP, int nIND){

	int i, j;
	init_random();
	*sigma2_F = calloc(K, sizeof(double));
	for (i=0; i<K; i++) *(*(sigma2_F) + i) = rand_double(.1, 1);

        *sigma2_Lambda = calloc(K, sizeof(double));
	*rowSd = calloc(nSNP, sizeof(double));
        for (i=0; i<K; i++) *(*(sigma2_Lambda) + i) = rand_double(.1, 1);

	*Factors = calloc(nSNP*K, sizeof(double));
	for (i=0; i<K; i++){
		for (j=0; j<nSNP; j++){
			*(*(Factors) + i*nSNP + j) = rand_normal(0, *(*(sigma2_F) + i));
		}
	}
	tr(*Factors, K, nSNP);

	*Lambda = calloc(K*nIND, sizeof(double));
	for (i=0; i<K*nIND; i++) *(*(Lambda) + i) = rand_normal_r();

	*Piz = calloc(K + 1, sizeof(double));
	*(*(Piz) + 0) = 1 - alpha;
	for (i=1; i<(K + 1); i++) *(*(Piz) + i) = alpha/K;

	*Pz = calloc(nSNP*(K + 1), sizeof(double));
	*C2 = calloc(K, sizeof(double)*K);
	for (i=0; i<K; i++) (*(C2))[i] = 2;

	return 0;
}

double smartInit(double *Lambda, double *sigma2_F, int nIND, int K, char *subSampleName){
	double *FactorsInit;
	double *GenoInit;
	int *missingInit;
	double *AF;
	double *MyCov;
	double err;
	int nrow, ncol, i, j, k, tmp;
	DataCounter(subSampleName, &nrow, &ncol, "genotypes");
	GenoInit = calloc(nrow*ncol, sizeof(double));
	if (extractData(subSampleName, nrow, ncol, GenoInit, "genotypes")) {
		printf("Invalid Initialisation File: %s with %d SNPs\n", subSampleName, ncol);
		free(GenoInit);
		return .0;
	} else {
		if (nIND == nrow){
			tr(GenoInit, nrow, ncol);
                	tmp = ncol;
        	        ncol = nrow;
	                nrow = tmp;
		}
        	FactorsInit = calloc(K*nrow, sizeof(double));
	        AF = calloc(nrow, sizeof(double));
		MyCov = calloc(nIND*nIND, sizeof(double));
		printf("Initialisation with data in file %s\n", subSampleName);
		DataImputation(GenoInit, &missingInit, nrow, ncol, "naive");
	        rowMeans(GenoInit, nrow, ncol, AF);
	        for (i=0; i<nrow; i++){
        	        for (j=0; j<ncol; j++) GenoInit[i*nIND + j] -= AF[i];
        	}
        	free(AF);
		covariance(GenoInit, MyCov, nrow, ncol, 1);	
		diagonalize(MyCov, ncol, K, sigma2_F, Lambda);
		tr(Lambda, nIND, K);

		free(GenoInit);
		free(MyCov);
		free(missingInit);
		free(FactorsInit);
		return err;
	}
}

int runMCMC(double *G, int *missing, double *Factors, double *sigma2_F, double *Lambda, double *Pz, double *Piz, double *sigma2, int K, double alpha, double Tau2, double *C2, double *sigma2_Lambda, double Isingb, int nSNP, int nIND, int nsteps, int burnin, int na, int prior, int runSVD, char *OutputFileName, char *subSampleName, double prop){

	int step, noutlier, i, j, k, no, updates;
	int *Z = calloc(nSNP, sizeof(int));
	double *Pz_chain = calloc(nSNP*(K + 1), sizeof(double));
	double *FLambda = calloc(nIND*nSNP, sizeof(double));
 	double *Lambda_chain = calloc(nIND*K, sizeof(double));
 	double *Factors_chain = calloc(nSNP*K, sizeof(double));
	double sqerr, sqerr_chain = 0, alpha_chain = 0;
	double *sigma2_F_chain = calloc(K, sizeof(double));
	double *Tau2_chain = calloc(K, sizeof(double));
	double *AF = calloc(nSNP, sizeof(double));
	double *PriorOdds = calloc(nSNP, sizeof(double));
        double *PO = calloc(nSNP, sizeof(double));

        rowMeans(G, nSNP, nIND, AF);
        for (i=0; i<nSNP; i++){
                for (j=0; j<nIND; j++) G[i*nIND + j] -= AF[i];
        }
	free(AF);
	if (runSVD){
	        sqerr = smartInit(Lambda, sigma2_F, nIND, K, subSampleName);
	        printf("Singular Value Decomposition done\n");
	}

	for (i=0; i<K; i++) sigma2_Lambda[i] = 1.0;

	for (step=0; step<nsteps; step++){

                updateFactors(Factors, Z, Lambda, G, Tau2, C2, sigma2_F, *sigma2, nSNP, nIND, K, prior);
                updates = updateSigma2_F(sigma2_F, Factors, *sigma2, Z, Tau2, C2, K, nSNP, prior);
                updateSigma2(sigma2, G, Factors, Lambda, sigma2_F, C2, Z, nSNP, nIND, K, prior);
		updatePz(Pz, Factors, alpha, sigma2_F, Z, Tau2, C2, *sigma2, Isingb, nSNP, K, prior, PriorOdds, PO, step, burnin);
		noutlier = 0;
		for (i=0; i<nSNP; i++) if (Z[i]) noutlier++;
		updateAlpha(&alpha, noutlier, nSNP);
		updateTau2(&Tau2, C2, sigma2_F, Factors, *sigma2, Z, K, nSNP, prior);
                updateLambda(Lambda, G, Factors, *sigma2, sigma2_Lambda, K, nSNP, nIND);
//		updateSigma2_Lambda(sigma2_Lambda, Lambda, K, nIND);

		if (step >= burnin){
			for (i=0; i<nSNP*(K + 1); i++) Pz_chain[i] += Pz[i]/(nsteps - burnin);
			for (i=0; i<nSNP*K; i++) Factors_chain[i] += Factors[i]/(nsteps - burnin);
			for (i=0; i<nIND*K; i++) Lambda_chain[i] += Lambda[i]/(nsteps - burnin);
			sqerr = 0;
			prodMatrix(Factors, Lambda, FLambda, nSNP, K, K, nIND);
			for (i=0; i<nSNP*nIND; i++) sqerr += (G[i] - FLambda[i])*(G[i] - FLambda[i])/(nIND*nSNP);
			sqerr_chain += sqerr/(nsteps - burnin);
			alpha_chain += alpha/(nsteps - burnin);
                        for (i=0; i<K; i++) Tau2_chain[i] += C2[i]/(nsteps - burnin);
                        for (i=0; i<K; i++) sigma2_F_chain[i] += sigma2_F[i]/(nsteps - burnin);
			if (!((step + 1)%10)){
				printf("MCMC step: %i\n", step + 1);
				printf("\tmean squarred error: %f\n", sqerr);
				printf("\talpha value: %f\n\tsigma2: %f\n\toutliers: %i\n\tsigma_F: ", alpha, *sigma2, noutlier);
				displayMatrix(sigma2_F, 1, K);
				printf("\tC2 values: ");
				displayMatrix(C2, 1, K);
				printf("\n");
			}

		}
		int nan = 0;
		for (i=0; i<K; i++) nan += isnan(sigma2_F[i]);
		nan += isnan(*sigma2);
		if (nan){step = nsteps; printf("Not a number in the chain, STOP\nPlease rerun PCAdapt\n");} 
	}

	double ConstantPriorO = ConstantPriorOdds();

	for (i=0; i<nSNP; i++){
//		PriorOdds[i] /= (nsteps - burnin);
		PriorOdds[i] = ConstantPriorO;
//		PriorOdds[i] = log10(PriorOdds[i]);
		
		PO[i] /= (nsteps - burnin);
//		PO[i] = log10(PO[i]);
/* Average Bayes Factor */
		PriorOdds[i] = PO[i] - PriorOdds[i];
	}
	int *order_v = malloc(sizeof(int)*K);
	order(sigma2_F_chain, order_v, K);
	writeResults(OutputFileName, Pz_chain, PriorOdds, PO, nSNP, K, sqerr_chain, alpha_chain, Tau2_chain, sigma2_F_chain, order_v);
	writeScores(OutputFileName, Lambda_chain, Factors_chain, K, nIND, nSNP, order_v);
	writeBF(OutputFileName, PriorOdds, Pz_chain, nSNP, K, prop, order_v);

	free(Z);
	free(PO);
	free(PriorOdds);
	free(FLambda);
	free(Pz_chain);
	free(Lambda_chain);
	free(Factors_chain);
	return 0;

}
