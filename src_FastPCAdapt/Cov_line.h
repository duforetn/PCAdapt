#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../src/matrix.h"
#include "../src/Data.h"
#include "Data__f.h"
#include "../src/linAlgebra.h"

#define NA 9
#define TOL .0001

/* get_row
 * Read a block of lines (snp), scale it and store it.
 */
int get_row(double *Geno, FILE *GenoFile, int nIND, double *mean, double *SNPSd, int sc, int blocksize);

/* add_to_cov
 * Given a block of snps, computes correlation matrix of the block in scratchCov (allocated) and add it to the existing matrix Cov.
 */
void add_to_cov(double *Cov, double *scratchCov, int nIND, double *Geno, int blocksize);

/* Cov_line
 * function to learn covariance (or cor) matrix in COv, while reading a data file in GenoFileName, by storing blocks of lines. The function also handles the case where several input files are used.
 *TODO: parallelize
 */
int Cov_line(double *Cov, double *SNPSd, int nSNP, int *nSNP_file, int nIND, int sc, char **GenoFileName, int nfile);
