/*
 * B = (X^T M^-1 X)^-1 X^T M^-1 Y
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//#include "utils.h"
#include "options.h"
#include "lapack.h"
#include "blas.h"
//#include "bio.h"

void bio_chol(int m, int n, int p, int t, 
              double *B, double *X, double *Phi, double *y, 
              double *h)
{
    double *xTSx,
           *M,
           *copyX, *copyY,
           ONE = 1.0,
           ZERO = 0,
           hh; 
    int    info,
           mp = m*p,
           nn = n*n,
           iONE = 1,
           i, j;

    double *tmp;

    xTSx  = (double *) malloc ( p * p * sizeof(double) );
    M     = (double *) malloc ( n * n * sizeof(double) );
    copyX = (double *) malloc ( m * p * n * sizeof(double) );
    copyY = (double *) malloc ( n * sizeof(double) );
    if ( xTSx == NULL || M == NULL || copyX == NULL )
    {
        fprintf(stderr, "[bio_var3] Not enough memory\n");
        exit(-1);
    }

    

    memcpy( copyX, X, m * p * n * sizeof(double) );
    for ( j = 0; j < t; j++ )
    {
        hh = h[j]*h[j];
        memcpy( copyX, X, m * p * n * sizeof(double) );
        memcpy( copyY, &y[j*n], n * sizeof(double) );
        memcpy( M, Phi, n * n * sizeof(double) );
        /* 1) M := h^2 Phi + (1-h^2) I) */
        dscal_(&nn, &hh, M, &iONE);
        for ( i = 0; i < n; i++ )
            M[i*n + i] += 1 - hh;

        /* 2) L * L' = M */
        dpotrf_(LOWER, &n, M, &n, &info);

        /* 3) X := inv(L) X */
        dtrsm_(LEFT, LOWER, NO_TRANS, NON_UNIT, &n, &mp, &ONE, M, &n, copyX, &n);

        /* 4) y := inv(L) y */
        dtrsv_(LOWER, NO_TRANS, NON_UNIT, &n, M, &n, copyY, &iONE);

        /* 5) B := xTSy := X' * y */
        dgemv_(TRANS, &n, &mp, &ONE, copyX, &n, copyY, &iONE, &ZERO, &B[j*mp], &iONE);

        for ( i = 0; i < m; i++ )
        {
            /* 6) xTSx := Xi' * Xi */
            dsyrk_(LOWER, TRANS, &p, &n, &ONE, &copyX[n*i*p], &n, &ZERO, xTSx, &p);

            /* 7) B = inv(xTSx) * xTSy */
            dposv_(LOWER, &p, &iONE, xTSx, &p, &B[j*mp + i*p], &p, &info);
            if (info != 0)
            {
                fprintf(stderr, "Error executing dposv: %d\n", info);
                exit(-1);
            }
        }
    }

    free(xTSx);
    free(M);
    free(copyX);
    free(copyY);
}    
