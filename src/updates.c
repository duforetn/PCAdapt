/*
   PCAdapt updates.c
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
#include "matrix.h"
#include "random.h"
#include "updates.h"
#include "linAlgebra.h"

#define TAU2INF 2.0
#define TAU2SUP 10.0
#define LOGALPHAINF -4.0
#define LOGALPHASUP -2.0

void updateMu(double *mu, double *G, double *Factors, double *Lambda, double sigma2, int nSNP, int nIND, int K){

	int i, j;
	double mu_i;
        double *meanMu = malloc(sizeof(double)*nSNP);
        double *diff = malloc(sizeof(double)*nIND*nSNP); 
	double *prod = malloc(sizeof(double)*nIND*nSNP);
        prodMatrix(Factors, Lambda, prod, nSNP, K, K, nIND);
        difMatrix(G, prod, diff, nSNP, nIND);
        rowMeans(diff, nSNP, nIND, meanMu);
        for (i=0; i<nSNP; i++){
		mu_i = rand_normal(meanMu[i], sigma2/nIND);
		for (j=0; j<nIND; j++) mu[i*nIND + j] = mu_i;
	}
	free(prod);
        free(meanMu);
        free(diff);

}

void updateLambda(double *Lambda, double *G, double *Factors, double sigma2, double *sigma2_Lambda, int K, int nSNP, int nIND){

	int i, j;
	double *meanLambda_j = calloc(K, sizeof(double));
	double *invSigmaLambda = calloc(K*K, sizeof(double));
	double *CholinvSigmaLambda = calloc(K*K, sizeof(double));
	double *trFactors = calloc(K*nSNP, sizeof(double));
	double *Gj = calloc(nSNP, sizeof(double));
	double *FGj = calloc(K, sizeof(double));
	double *rLambda_j = calloc(K, sizeof(double));
	double *trLambda = calloc(K*nIND, sizeof(double));
	
	for (i=0; i<(nSNP*K); i++) trFactors[i] = Factors[i];
	tr(trFactors, nSNP, K);

	prodMatrix(trFactors, Factors, invSigmaLambda, K, nSNP, nSNP, K);
	for (i=0; i<K; i++){
		for (j=0; j<K; j++){
			invSigmaLambda[i*K + j] = invSigmaLambda[i*K + j]/sigma2;
		}
		invSigmaLambda[i*K + i] += 1.0/sigma2_Lambda[i];
	}
	invSPDMatrix(invSigmaLambda, K);
	cholesky(invSigmaLambda, K, CholinvSigmaLambda, 1);

	for (j=0; j<nIND; j++){
		for (i=0; i<nSNP; i++) Gj[i] = G[i*nIND + j]/sigma2;
		prodMatrix(trFactors, Gj, FGj, K, nSNP, nSNP, 1);
		prodMatrix(invSigmaLambda, FGj, meanLambda_j, K, K, K, 1);

		mvn_rand(meanLambda_j, CholinvSigmaLambda, K, rLambda_j);
		for (i=0; i<K; i++) trLambda[i*nIND + j] = rLambda_j[i];

	}
	for (i=0; i<K*nIND; i++) Lambda[i] = trLambda[i];
	//tr(Lambda, nIND, K);

	free(trLambda);
	free(meanLambda_j);
	free(invSigmaLambda);
	free(CholinvSigmaLambda);
	free(trFactors);
	free(Gj);
	free(FGj);
	free(rLambda_j);
}
double ConstantPriorOdds(){
	int i;
        double PriorDraws[100];
        for (i=0; i<100; i++) PriorDraws[i] = pow(10, rand_double(LOGALPHAINF, LOGALPHASUP)); 
        double ConstantPriorO = 0;
        for (i=0; i<100; i++) ConstantPriorO += (log10(PriorDraws[i]) - log10(1 - PriorDraws[i]))/100; 
	return ConstantPriorO;

}

void updatePz(double *Pz, double *Factors, double alpha, double *sigma2_F, int *Z, double Tau2, double *C2, double sigma2, double Isingb, int nSNP, int K, int prior, double *PriorOdds, double *PO, int step, int burnin){

	int i, j, k, kmax = 0;
	double c, L;
	double sumPiz = 0;
	double Pkmax = 0, Lmax = 0;
	double *Pz_i = calloc(K + 1, sizeof(double));
	double Piz[K + 1];
	
	if (prior == 0) sigma2 = 1.0;
	Piz[0] = 1 - alpha;
	for (j=1; j<K + 1; j++) Piz[j] = alpha/K;

	for (i=0; i<nSNP; i++){
		Piz[0] = 1 - alpha;
		for (j=1; j<K + 1; j++) Piz[j] = alpha/K;
		if(i>0) Piz[Z[i-1]] *= exp(Isingb);
		if(i<(nSNP - 1)) Piz[Z[i+1]] *= exp(Isingb);
		
		Pz[i*(K + 1)] = Piz[0];
		c = Piz[0];
		Pkmax = 0;
		Lmax = 0;
		kmax = 1;
		for (k=1; k<K + 1; k++){
			Pz[i*(K + 1) + k] = Piz[k];
//CHANGE 1
			L = exp(Factors[i*K + k - 1]*Factors[i*K + k - 1]*(C2[k - 1] - 1.0)/(C2[k - 1]*2.0*sigma2*(sigma2_F[k - 1])))/sqrt(C2[k - 1]);
			Pz[i*(K + 1) + k] *= L/Piz[0];
			c += Pz[i*(K + 1) + k];
                        if(Pkmax < Pz[i*(K + 1) + k]) kmax = k;
			if(Pkmax < Pz[i*(K + 1) + k]) Pkmax = Pz[i*(K + 1) + k];
		}
//		Pkmax = Pz[i*(K + 1) + kmax];
//		Lmax = Pz[i*(K + 1) + kmax]/Piz[kmax];
                if (step >= burnin) PriorOdds[i] += log10(Piz[kmax]/Piz[0]);
//		sumPiz = 0;
//		for (k=1; k<K + 1; k++) sumPiz += Piz[k];
		if (step >= burnin) PO[i] += log10(Pz[i*(K + 1) + kmax]);

		if (!isnan(c) && !isinf(c)){
			for (k=0; k<K + 1; k++) Pz[i*(K + 1) + k] = Pz[i*(K + 1) + k]/c;
		} else {
			/* Case of extreme values */
			int imax = 0;
			for (k=0; k<K; k++) if(Factors[i*K + k]*Factors[i*K + k]/sigma2_F[k] > Factors[i*K + kmax]*Factors[i*K + kmax]/sigma2_F[kmax]) imax = k;
			for (k=0; k<K + 1; k++) Pz[i*(K + 1) + k] = 0;
			Pz[i*(K + 1) + imax + 1] = 1;
		}
                for (k=0; k<K + 1; k++) Pz_i[k] = Pz[i*(K + 1) + k];
                Z[i] = rand_vector(Pz_i, K + 1);
	}
	free(Pz_i);
}

void updateAlpha(double *alpha, int noutlier, int nSNP){

	double v = (double)((LOGALPHASUP-LOGALPHAINF)/8.0)*((LOGALPHASUP-LOGALPHAINF)/8.0);
	double prop_alpha = *alpha*pow(10.0, rand_normal(.0, v));
	if (!(log10(prop_alpha) > LOGALPHASUP || log10(prop_alpha) < LOGALPHAINF)){
		double ratio =((double) noutlier)*(log(prop_alpha) - log(*alpha)) + ((double) (nSNP - noutlier))*(log(1.0 - prop_alpha) - log(1.0 - *alpha));
		double tg = ((double) noutlier)*(log(prop_alpha) - log(*alpha));
		double td = ((double) (nSNP - noutlier))*(log(1.0 - prop_alpha) - log(1.0 - *alpha));
	//printf("alpha: %f, Proposal: %f, no:%i, logratio: (%f) (%f)\n", *alpha, prop_alpha, noutlier, tg, td);
		if (log(drand()) < ratio) *alpha = prop_alpha;
	}
}

void updateFactors(double *Factors, int *Z, double *Lambda, double *G, double Tau2, double *C2, double *sigma2_F, double sigma2, int nSNP, int nIND, int K, int prior){

	int i, j, l, ii;
	double *invSigmaFactors = calloc(K*K, sizeof(double));
	double *trLambda = calloc(K*nIND, sizeof(double));
	double *LtrL = calloc(K*K, sizeof(double));
	double *CholinvSigmaFactors = calloc(K*K, sizeof(double));
        double *invSigmaFactors_outlier = calloc(K*K, sizeof(double));
        double *CholinvSigmaFactors_outlier = calloc(K*K, sizeof(double));
	double *Gi = calloc(nIND, sizeof(double));
	double *LGi = calloc(K, sizeof(double));
	double *meanFactor_i = calloc(K, sizeof(double));
        double *rFactor_i = calloc(K, sizeof(double));

        for (i=0; i<(K*nIND); i++) trLambda[i] = Lambda[i];
        tr(trLambda, K, nIND);
	prodMatrix(Lambda, trLambda, LtrL, K, nIND, nIND, K);
	for (ii=0; ii<K; ii++){
		for (j=0; j<K; j++){
//			invSigmaFactors[i*K + j] = LtrL[i*K + j]/sigma2;
			invSigmaFactors[ii*K + j] = LtrL[ii*K + j];
		}
		invSigmaFactors[ii*K + ii] += 1.0/sigma2_F[ii];
	}
	invSPDMatrix(invSigmaFactors, K);
	cholesky(invSigmaFactors, K, CholinvSigmaFactors, 1);
        for (i=0; i<nSNP; i++){

		if (Z[i] == 0){
//			for (j=0; j<nIND; j++) Gi[j] = G[i*nIND + j]/sigma2;
			for (j=0; j<nIND; j++) Gi[j] = G[i*nIND + j];
			prodMatrix(Lambda, Gi, LGi, K, nIND, nIND, 1);
			prodMatrix(invSigmaFactors, LGi, meanFactor_i, K, K, K, 1);
			mvn_rand(meanFactor_i, CholinvSigmaFactors, K, rFactor_i);
			for (j=0; j<K; j++) Factors[i*K + j] = rFactor_i[j];
		} else {
			/* change the posterior variance matrix with Variance inflation */
		        for (l=0; l<K; l++){
		                for (j=0; j<K; j++){
//                	        	invSigmaFactors_outlier[l*K + j] = LtrL[l*K + j]/sigma2;
					invSigmaFactors_outlier[l*K + j] = LtrL[l*K + j];
                		}
				if (l == (Z[i] - 1)){
					invSigmaFactors_outlier[l*K + l] += 1.0/(sigma2_F[l]*C2[Z[i] - 1]);
				} else {
					invSigmaFactors_outlier[l*K + l] += 1.0/(sigma2_F[l]);
				}
		        }
			invSPDMatrix(invSigmaFactors_outlier, K);
		        cholesky(invSigmaFactors_outlier, K, CholinvSigmaFactors_outlier, 1);
//CHANGE 2
			if (prior == 0) for (j=0; j<nIND; j++) Gi[j] = G[i*nIND + j]/(sigma2)*sigma2;
			for (j=0; j<nIND; j++) Gi[j] = G[i*nIND + j];
                        prodMatrix(Lambda, Gi, LGi, K, nIND, nIND, 1);
                        prodMatrix(invSigmaFactors_outlier, LGi, meanFactor_i, K, K, K, 1);
                        mvn_rand(meanFactor_i, CholinvSigmaFactors_outlier, K, rFactor_i);
                        for (j=0; j<K; j++){
				Factors[i*K + j] = rFactor_i[j];
			}
		}
	}

	free(Gi);
	free(LGi);
	free(LtrL);
	free(trLambda);
	free(invSigmaFactors);
        free(CholinvSigmaFactors);
        free(invSigmaFactors_outlier);
        free(CholinvSigmaFactors_outlier);
	free(meanFactor_i);
	free(rFactor_i);
}

int updateSigma2_F(double *sigma2_F, double *Factors, double sigma2, int *Z, double Tau2, double *C2, int K, int nSNP, int prior){

	int i, k;
	double s;
	int updates = K;
	for (k=0; k<K; k++){
		s = 0;
		for (i=0; i<nSNP; i++){
			if (Z[i] == (k + 1)){
				s += (Factors[i*K + k])*(Factors[i*K + k])/C2[Z[i] - 1];
			} else {
				s += (Factors[i*K + k])*(Factors[i*K + k]);
			}
		}
//CHANGE 3
		if (prior == 1) s /= sigma2;
			if (s > .00001/nSNP && s < 100*nSNP)
		//		sigma2_F[k] = rand_invgamma(nSNP/2.0, s/2.0);
				sigma2_F[k] = s/rand_chisq(nSNP);
			else updates--;
	}
	return updates;
}

void updateSigma2_Lambda(double *sigma2_Lambda, double *Lambda, int K, int nIND){

        int i, k;
        double s;
        int updates = K;
        for (k=0; k<K; k++){
                s = 0;
                for (i=0; i<nIND; i++){
                        s += (Lambda[k*nIND + i])*(Lambda[k*nIND + i]);
                }
//                if (s > .1/nIND && s < 100*nIND)
//                        sigma2_Lambda[k] = s/rand_chisq(nIND);
		sigma2_Lambda[k] = s/nIND;
        }
}

void updateSigma2(double *sigma2, double *G, double *Factors, double *Lambda, double *sigma2_F, double *C2, int *Z, int nSNP, int nIND, int K, int prior){
	
	int i, j;
	double s = 0;
	double *FLambda = calloc(nSNP*nIND, sizeof(double));
	prodMatrix(Factors, Lambda, FLambda, nSNP, K, K, nIND);

	for (i=0; i<nSNP; i++){
		for (j=0; j<nIND; j++){
			s += (G[i*nIND + j] - FLambda[i*nIND + j])*(G[i*nIND + j] - FLambda[i*nIND + j]);
		}
	}
//CHANGE 0
	if (prior == 1){
        	for (i=0; i<nSNP; i++){
	                for (j=0; j<K; j++){
                	        if (Z[i] == (j + 1)){
                        	        s += (Factors[i*K + j])*(Factors[i*K + j])/(sigma2_F[j]*C2[Z[i] - 1]);
	                        } else {
        	                        s += (Factors[i*K + j])*(Factors[i*K + j])/sigma2_F[j];
				}
                        }
	        }
	}
//	*sigma2 = rand_invgamma(nSNP*nIND/2, s/2.0);
	*sigma2 = s/rand_chisq(nSNP*nIND + K*nSNP);
	free(FLambda);
}

void updateTau2(double *Tau2, double *C2, double *sigma2_F, double *Factors, double sigma2, int *Z, int K, int nSNP, int prior){

	int i, k, f;
        double Tau2_prop;
	for (f=0; f<K; f++){
		Tau2_prop = rand_normal(C2[f], 1.0);
		if (!(Tau2_prop > TAU2SUP || Tau2_prop < TAU2INF)){
			double p_old = 0, p_prop = 0;
			if (prior == 0) sigma2 = 1.0;

			for (i=0; i<nSNP; i++){
				if(Z[i] == f + 1){
					k = Z[i] - 1;
//CHANGE 4
					p_old += (-Factors[i*K + k]*Factors[i*K + k]/(2.0*(C2[f])*sigma2_F[k]*sigma2)) - log(C2[f])/2.0;
        		                p_prop += (-Factors[i*K + k]*Factors[i*K + k]/(2.0*Tau2_prop*sigma2_F[k]*sigma2)) - log(Tau2_prop)/2.0;
				}
			}
//	printf("Tau2 :%f prop: %f, ratio: %f = %f - %f\n", *Tau2, Tau2_prop, p_prop - p_old, p_prop, p_old);
			if(log(drand()) < (p_prop - p_old)) C2[f] = Tau2_prop;
		}
	}
}
