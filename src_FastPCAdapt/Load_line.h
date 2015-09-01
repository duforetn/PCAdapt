#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../src/matrix.h"
#include "../src/Data.h"
#include "../src/linAlgebra.h"
#include "Data__f.h"

#define NA 9

//Load_line
//computes the U matrix, by block of snps, with no storage.

int Load_line(double *U, double *Sigma, double *V, double *miss, double *mAF, int nSNP, int *nSNP_file, int K, int nIND, int sc, char **GenoFileName, int nfile);
