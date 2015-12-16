#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../src/matrix.h"
#include "../src/Data.h"
#include "Data__f.h"
#include "../src/linAlgebra.h"

#define NA 9
#define TOL .0001

//get_row
//read a block of lines (snp), scale it and store it


int get_row(double *Geno, FILE *GenoFile, int nIND, double *mean, double *SNPSd, int *pairwiseObs, int sc, int blocksize, int calculateCov);

//add_to_cov
//Given a block of snps, computes correlation matrix of the block and add it to the existing matrix.

void add_to_cov(double *Cov, double *scratchCov, int nIND, double *Geno, int blocksize);

//add_to_obs
void add_to_obs(int *pairwiseObs, double *scratchCov, int nIND, double *Geno, int blocksize);

//Cov_line
//function to learn covariance (or cor) matrix, while reading a data file, by storing blocks of lines.
//TODO: parallelize

int Cov_line(double *Cov, double *SNPSd, int nSNP, int *nSNP_file, int nIND, int *pairwiseObs, int sc, char **GenoFileName, int nfile);
