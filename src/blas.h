#ifndef __BLAS_H
#define __BLAS_H

/*
 * BLAS 1
 */

extern void   drotg_ ( double *a,  double *b,  double *c, double *s );
extern void   drotmg_( double *d1, double *d2, double *a, double *b, double *param );
extern void   drot_  ( int *n, double *x, int *incx, double *y, int *incy, double *c, double *s );
extern void   drotm_ ( int *n, double *x, int *incx, double *y, int *incy, double *param );
extern void   dswap_ ( int *n,                double *x, int *incx, double *y, int *incy );
extern void   dscal_ ( int *n, double *alpha, double *x, int *incx );
extern void   dcopy_ ( int *n,                double *x, int *incx, double *y, int *incy );
extern void   daxpy_ ( int *n, double *alpha, double *x, int *incx, double *y, int *incy );
extern double ddot_  ( int *n,                double *x, int *incx, double *y, int *incy );
extern double ddotu_ ( int *n,                double *x, int *incx, double *y, int *incy );
extern double ddotc_ ( int *n,                double *x, int *incx, double *y, int *incy );
extern double dnrm2_ ( int *n,                double *x, int *incx                       );
extern double dasum_ ( int *n,                double *x, int *incx                       );
extern int    idamax_( int *n,                double *x, int *incx                       );

/*
 * BLAS 2
 */

extern void dgemv_( char *trans, 
					int *m, int *n, 
					double *alpha, double *A, int *ldA, double *x, int *incx, 
					double *beta,                       double *y, int *incy );
extern void dgbmv_( char *trans,
					int *m, int *n, int *kl, int *ku, 
					double *alpha, double *A, int *ldA, double *x, int *incx, 
					double *beta,                       double *y, int *incy );
extern void dsymv_( char *uplo,
					int *n,
					double *alpha, double *A, int *ldA, double *x, int *incx, 
					double *beta,                       double *y, int *incy );
extern void dsbmv_( char *uplo,
					int *n, int *k,
					double *alpha, double *A, int *ldA, double *x, int *incx, 
					double *beta,                       double *y, int *incy );
extern void dhpmv_( char *uplo,
					int *n,
					double *alpha, double *AP, double *x, int *incx,
					double *beta,              double *y, int *incy );
extern void dtrmv_( char *uplo, char *trans, char *diag,
					int *n,
					double *A, int *ldA, double *x, int *incx );
extern void dtbmv_( char *uplo, char *trans, char *diag,
					int *n, int *k,
					double *A, int *ldA, double *x, int *incx );
extern void dtpmv_( char *uplo, char *trans, char *diag,
					int *n,
					double *AP, double *x, int *incx );
extern void dtrsv_( char *uplo, char *trans, char *diag,
					int *n,
					double *A, int *ldA, double *x, int *incx );
extern void dtbsv_( char *uplo, char *trans, char *diag,
					int *n, int *k,
					double *A, int *ldA, double *x, int *incx );
extern void dtpsv_( char *uplo, char *trans, char *diag,
					int *n,
					double *AP, double *x, int *incx );

extern void dger_ ( int *m, int *n,
					double *alpha, double *x, int *incx, 
								   double *y, int *incy, 
					double *A, int *ldA );
extern void dsyr_ ( char *uplo,
					int *n, 
					double *alpha, double *x, int *incx,
					double *A, int *ldA );
extern void dspr_ ( char *uplo,
					int *n, 
					double *alpha, double *x, int *incx,
					double *AP );
extern void dsyr2_( char *uplo,
					int *n, 
					double *alpha, double *x, int *incx,
								   double *y, int *incy,
					double *A, int *ldA );
extern void dspr2_( char *uplo,
					int *n, 
					double *alpha, double *x, int *incx,
								   double *y, int *incy,
					double *AP );

/*
 * BLAS 3
 */

extern void dgemm_( char *transA, char *transB, 
		            int *m, int *n, int *k, 
					double *alpha, double *A, int *ldA, double *B, int *ldB, 
					double *beta,  double *C, int *ldC );
extern void dsymm_( char *side, char *uplo, 
		            int *m, int *n, 
					double *alpha, double *A, int *ldA, double *B, int *ldB, 
					double *beta,  double *C, int *ldC );
extern void dsyrk_( char *uplo, char *trans, 
		            int *n, int *k, 
					double *alpha, double *A, int *ldA,
					double *beta,  double *C, int *ldC );
extern void dsyr2k_( char *uplo, char *trans, 
		             int *n, int *k, 
					 double *alpha, double *A, int *ldA, double *B, int *ldB, 
					 double *beta,  double *C, int *ldC );
extern void dtrmm_(char *side, char *uplo, char *transA, char *diag, 
				   int *m, int *n, 
				   double *alpha, double *A, int *ldA, 
				   double *B, int *ldB);
extern void dtrsm_(char *side, char *uplo, char *transA, char *diag, 
				   int *m, int *n, 
				   double *alpha, double *A, int *ldA, 
				   double *B, int *ldB);

#endif // __BLAS_H
