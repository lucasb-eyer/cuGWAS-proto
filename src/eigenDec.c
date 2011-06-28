#include "blas.h"
#include "lapack.h"
#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
void eigenDec(int n, double *Phi, double *Z, double *W)
{
    int nb = 192;
    int idummy, nCompPairs, *isuppz, *iwork, info,
        lwork = n * (nb + 6),
        liwork = 10 * n;
    double ddummy = -1.0, *work;

    work = (double *) malloc ( lwork * sizeof(double) );
    iwork = (int *) malloc ( liwork * sizeof(int) );
    isuppz = (int *) malloc ( 2 * n * sizeof(int) );

    dsyevr_("V", "A", "L", &n, Phi, &n, &ddummy, &ddummy, &idummy, &idummy, 
    &ddummy, &nCompPairs, W, Z, &n, isuppz, work, &lwork, iwork, &liwork,
    &info);

    if (info != 0)
    {
        fprintf(stderr, "Error executing dsyevr\n");
        exit(-1);
    }
}

