/*
 * B = (X^T M^-1 X)^-1 X^T M^-1 Y
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

//#include "utils.h"
//#include "bio.h"
#include "options.h"
#include "lapack.h"
#include "blas.h"

/* 
 * Z eigenvectors of Phi
 * W eigenvalues of Phi
 */

/* we assume p << n reusing work in dsysv */
void bio_eigen(int m, int n, int p, int t, 
               double *B, double *X, double *Z, double *W, double *y, double *h)
{
    double *ZtX,     // Z' X
           *ZtY,     // Z' Y
           *Winv,    // inv( h^2 W - (1 - h^2) I )
           *XtZWinv, // X' Z inv(sqrt(W))
           *xtSx,    // x_i' inv(S) x_i
           ONE = 1.0,
           ZERO = 0,
           nrm;
    int    info,
           mp = m*p,
           iONE = 1,
           i, j, k, l;

    ZtX = (double *) malloc ( mp * n * sizeof(double) );
    ZtY = (double *) malloc ( n * t * sizeof(double) );
    Winv = (double *) malloc ( n * sizeof(double) );
    XtZWinv = (double *) malloc ( mp * n * sizeof(double) );
    xtSx = (double *) malloc ( p * p * sizeof(double) );

    if (xtSx == NULL)
    {
        fprintf(stderr, "Not enough memory\n");
        exit(-1);
    }
    printf("It fits\n");
    /*printf("m: %4d\n", m);*/
    /*printf("n: %4d\n", n);*/
    /*printf("p: %4d\n", p);*/
    /*printf("t: %4d\n", t);*/

    /* 3) ZtX <- X^T * Z */
    dgemm_("T", "N", &mp, &n, &n, &ONE, X, &n, Z, &n, &ZERO, ZtX, &mp);

    /* 6) ZtY = Z' * Y */
    dgemm_("T", "N", &n, &t, &n, &ONE, Z, &n, y, &n, &ZERO, ZtY, &n);


    for ( i = 0; i < t; i++ )
    {
        /* 2) W = sqrt(alpha W - beta I)^-1 */
        // Best order? sqrt - inv
        for ( k = 0; k < n; k++)
            Winv[k] = sqrt(1.0 / (h[k]*h[k] * W[k] + (1 - h[k]*h[k])));

        /* X' * Z  * sqrt(Winv) */
        for (l = 0; l < n*mp; l++) XtZWinv[l] = 0.0;
        for (l = 0; l < n; l++ )
            daxpy_(&mp, &Winv[l], &ZtX[l*mp], &iONE, &XtZWinv[l*mp], &iONE);

        /* sqrt(Winv) * ZtY */
        for (l = 0; l < n; l++ )
            /*dscal_(&t, &Winv[l], &ZtY[l], &n);*/
            ZtY[i*n + l] *= Winv[l];

        /* 7) y = XtZWinv * y */
        dgemv_("N", &mp, &n, &ONE, XtZWinv, &mp, &ZtY[i*n], &iONE, &ZERO, &B[i*mp], &iONE);

        for ( j = 0; j < m; j++ )
        {
            /* 5) W = XtZWinv * K^T */
            dsyrk_("L", "N", &p, &n, &ONE, &XtZWinv[j*p], &mp, &ZERO, xtSx, &p);

            /* 8) W^-1 * y */
            dposv_("L", &p, &iONE, xtSx, &p, &B[i*mp + j*p], &p, &info);
            if (info != 0)
            {
                fprintf(stderr, "Error executing dposv: %d\n", info);
                exit(-1);
            }
        }
    }

    free(ZtX);
    free(ZtY);
    free(Winv);
    free(XtZWinv);
    free(xtSx);
}    
