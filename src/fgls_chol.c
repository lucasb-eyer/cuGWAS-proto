#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <sys/types.h>
#include <sys/stat.h>

#include <sys/time.h>
#include <time.h>

#include <pthread.h>
#include <semaphore.h>
#include <aio.h>

#include "blas.h"
#include "lapack.h"
#include "options.h"
#include "common.h"
#include "io.h"
#include "timing.h"
#include "fgls_chol.h"

#if VAMPIR
	#include "vt_user.h"
#endif

void* compute_thread_func(void* in);

/*
 * Cholesky-based solution of the Feasible Generalized Least-Squares problem
 */
int fgls_chol(int n, int p, int m, int t, int wXL, int wXR,
				 int x_b, int y_b, int num_threads,
				 char *Phi_path, char *h2_path, char *sigma2_path,
				 char *XL_path, char *XR_path, char *Y_path,
				 char *B_path, char *V_path)
{
	FGLS_config_t cf;

	FILE *Phi_fp, *h_fp, *sigma_fp,
		 *XL_fp, *XR_fp, *Y_fp,
		 *B_fp, *V_fp;

	double *Phi;
	double *M;
	double *h;
	double *sigma;
	double alpha;
	double beta;

	double *xl, *xl_b, *xltxl;

	double *XL[2]; // XL and a copy (XL is overwritten at every iteration of j)
	double *XL_b;  // Top part of b ( in inv(S) b )
	double *XLtXL; // TopLeft part of S ( in inv(S) b )

	double *X[NUM_BUFFERS_PER_THREAD];
	double *Y[NUM_BUFFERS_PER_THREAD];
	double *B[NUM_BUFFERS_PER_THREAD];
	double *V[NUM_BUFFERS_PER_THREAD];

	double ONE = 1.0;
	double ZERO = 0.0;
	int iONE = 1;

	int ib, i, j, k;
	int nn = n * n;
	/*int x_inc, y_inc;*/
	double *Bij, *Vij;
	int info;

	/*printf("n: %d\np: %d\nm: %d\nt: %d\nwXL: %d\nwXR: %d\n", n, p, m, t, wXL, wXR);*/
	/*printf("x_b: %d\ny_b: %d\nnths: %d\n", x_b, y_b, num_threads);*/
	/*printf("Phi: %s\nh2: %s\ns2: %s\n", Phi_path, h2_path, sigma2_path);*/
	/*printf("XL: %s\nXR: %s\nY: %s\n", XL_path, XR_path, Y_path);*/
	/*printf("B: %s\nV: %s\n", B_path, V_path);*/

	initialize_config(
			&cf,
			n, p, m, t, wXL, wXR,
			x_b, y_b, num_threads,
			Phi_path, h2_path, sigma2_path,
			XL_path, XR_path, Y_path,
			B_path, V_path
	);

	Phi = ( double* ) fgls_malloc ( cf.n * cf.n * sizeof(double) );
	M   = ( double* ) fgls_malloc ( cf.n * cf.n * sizeof(double) );

	for (i = 0; i < NUM_BUFFERS_PER_THREAD; i++) 
	{
		X[i] = ( double* ) fgls_malloc ( cf.NUM_COMPUTE_THREADS * cf.x_b * cf.p * cf.n * sizeof(double) );
		Y[i] = ( double* ) fgls_malloc ( cf.NUM_COMPUTE_THREADS * cf.y_b * cf.n * sizeof(double) );
		B[i] = ( double* ) fgls_malloc ( cf.NUM_COMPUTE_THREADS * cf.x_b * cf.y_b * cf.p * sizeof(double) );
		V[i] = ( double* ) fgls_malloc ( cf.NUM_COMPUTE_THREADS * cf.x_b * cf.y_b * cf.p * cf.p * sizeof(double) );
	}
	// The original in XL[0]
	// The copy to overwrite in XL[1]
	XL[0] = ( double* ) fgls_malloc ( cf.wXL * cf.n * sizeof(double) );
	XL[1] = ( double* ) fgls_malloc ( cf.NUM_COMPUTE_THREADS * cf.wXL * cf.n * cf.y_b * sizeof(double) );
	XL_b  = ( double* ) fgls_malloc ( cf.NUM_COMPUTE_THREADS * cf.wXL * cf.y_b * sizeof(double) );
	XLtXL = ( double* ) fgls_malloc ( cf.NUM_COMPUTE_THREADS * cf.wXL * cf.wXL * cf.y_b * sizeof(double) );
	xl = XL[1];
	xl_b = XL_b;
	xltxl = XLtXL;

	XR_fp = fopen( cf.XR_path, "rb");
	Y_fp = fopen( cf.Y_path, "rb");
	B_fp = fopen( cf.B_path, "wb");
	V_fp = fopen( cf.V_path, "wb");

	h     = ( double* ) fgls_malloc ( cf.t * sizeof(double) );
	sigma = ( double* ) fgls_malloc ( cf.t * sizeof(double) );
	XL_fp = fopen( cf.XL_path, "rb" );
	sync_read( XL[0], XL_fp, cf.wXL * cf.n, 0 );
	fclose( XL_fp );
	h_fp = fopen( cf.h_path, "r");
	sync_read(h, h_fp, cf.t, 0);
	fclose( h_fp );
	sigma_fp = fopen( cf.sigma_path, "r");
	sync_read(sigma, sigma_fp, cf.t, 0);
	fclose( sigma_fp );
	
	/* Load Phi */
	Phi_fp = fopen( cf.Phi_path, "rb" );
	sync_read( Phi, Phi_fp, cf.n * cf.n, 0 );
	fclose( Phi_fp );

	double *x_cur	= X[0];
	double *x_next = X[1];
	double *y_cur	= Y[0];
	double *y_next = Y[1];
	double *b_cur	= B[0];
	double *b_prev = B[1];
	double *v_cur	= V[0];
	double *v_prev = V[1];

	struct aiocb aiocb_x_cur,	aiocb_x_next,
				 aiocb_y_cur,	aiocb_y_next,
				 aiocb_b_prev, aiocb_b_cur,
				 aiocb_v_prev, aiocb_v_cur;

	const struct aiocb ** aiocb_x_cur_l,//	= { &aiocb_x_cur }, 
		               ** aiocb_x_next_l,// = { &aiocb_x_next },
				       ** aiocb_y_cur_l,//	= { &aiocb_y_cur },
				       ** aiocb_y_next_l,// = { &aiocb_y_next },
				       ** aiocb_b_prev_l,// = {	aiocb_b_prev },
				       ** aiocb_b_cur_l,//	= {	aiocb_b_cur };
				       ** aiocb_v_prev_l,// = {	aiocb_v_prev },
				       ** aiocb_v_cur_l;//	= {	aiocb_v_cur };

	aiocb_x_cur_l  = (const struct aiocb **) fgls_malloc (sizeof(struct aiocb *));
	aiocb_x_next_l = (const struct aiocb **) fgls_malloc (sizeof(struct aiocb *));
	aiocb_y_cur_l  = (const struct aiocb **) fgls_malloc (sizeof(struct aiocb *));
	aiocb_y_next_l = (const struct aiocb **) fgls_malloc (sizeof(struct aiocb *));
	aiocb_b_prev_l = (const struct aiocb **) fgls_malloc (sizeof(struct aiocb *));
	aiocb_b_cur_l  = (const struct aiocb **) fgls_malloc (sizeof(struct aiocb *));
	aiocb_v_prev_l = (const struct aiocb **) fgls_malloc (sizeof(struct aiocb *));
	aiocb_v_cur_l  = (const struct aiocb **) fgls_malloc (sizeof(struct aiocb *));

	aiocb_x_cur_l[0]  = &aiocb_x_cur;
	aiocb_x_next_l[0] = &aiocb_x_next;
	aiocb_y_cur_l[0]  = &aiocb_y_cur;
	aiocb_y_next_l[0] = &aiocb_y_next;
	aiocb_b_prev_l[0] = &aiocb_b_prev;
	aiocb_b_cur_l[0]  = &aiocb_b_cur;
	aiocb_v_prev_l[0] = &aiocb_v_prev;
	aiocb_v_cur_l[0]  = &aiocb_v_cur;

	fgls_aio_read( &aiocb_x_cur,
			       fileno( XR_fp ), x_cur,
				   MIN( x_b, m ) * wXR * n * sizeof(double), 0 );

	fgls_aio_read( &aiocb_y_cur,
			       fileno( Y_fp ), y_cur,
				   n * sizeof(double), 0 );

	int iter = 0;
	for ( j = 0; j < t; j++ )
	{
		struct aiocb *aiocb_y_next_p = (struct aiocb *)aiocb_y_next_l[0];
		fgls_aio_read( aiocb_y_next_p,
					   fileno( Y_fp ), y_next,
					   j + 1 >= t ? 0 : n * sizeof(double),
					   j + 1 >= t ? 0 : (j+1) * n * sizeof(double) );

		memcpy( M, Phi, n * n * sizeof(double) );
		/* 1) M := sigma * ( h^2 Phi - (1 - h^2) I ) */
		alpha = h[j] * sigma[j];
		beta  = (1 - h[j]) * sigma[j];
		dscal_(&nn, &alpha, M, &iONE);
		for ( i = 0; i < n; i++ )
			M[i*n + i] = M[i*n + i] + beta;

		/* 2) L * L' = M */
		dpotrf_(LOWER, &n, M, &n, &info);
		if (info != 0)
		{
			char err[STR_BUFFER_SIZE];
			snprintf(err, STR_BUFFER_SIZE, "dpotrf(M) failed (info: %d)", info);
			error_msg(err, 1);
		}

		/* 3) WL := inv(L) XL */
		memcpy( xl, XL[0], wXL * n * sizeof(double) );
		dtrsm_(LEFT, LOWER, NO_TRANS, NON_UNIT, &n, &wXL, &ONE, M, &n, xl, &n);

		fgls_aio_suspend( aiocb_y_cur_l, 1, NULL );
		/* 4) y := inv(L) y */
		dtrsv_(LOWER, NO_TRANS, NON_UNIT, &n, M, &n, y_cur, &iONE);

		/* 5) B_t := XL' * y */
		dgemv_(TRANS, &n, &wXL, &ONE, xl, &n, y_cur, &iONE, &ZERO, xl_b, &iONE);
		/*for (int k = 1; k < m; k++)*/
		/*dcopy_(&widthXL, &B[j*mp], &iONE, &B[j*mp + k*p], &iONE);*/

		dsyrk_(LOWER, TRANS, &wXL, &n, &ONE, xl, &n, &ZERO, xltxl, &wXL);
		for (ib = 0; ib < m; ib += x_b) 
		{
			struct aiocb *aiocb_x_next_p = (struct aiocb *)aiocb_x_next_l[0];
			fgls_aio_read( aiocb_x_next_p,
					       fileno( XR_fp ), x_next,
						   (ib + x_b) >= m ? MIN( x_b, m ) * wXR * n * sizeof(double) : MIN( x_b * wXR * n, (m - (ib + x_b)) * wXR * n ) * sizeof(double),
						   (ib + x_b) >= m ? 0 : (ib + x_b) * wXR * n * sizeof(double) );

			fgls_aio_suspend( aiocb_x_cur_l, 1, NULL );
			/* 3) WR := inv(L) XR */
			int x_inc = MIN(x_b, m - ib);
			int rhss  = wXR * x_inc;
			dtrsm_(LEFT, LOWER, NO_TRANS, NON_UNIT, &n, &rhss, &ONE, M, &n, x_cur, &n);

			for (i = 0; i < x_inc; i++)
			{
				Bij = &b_cur[i * p];
				Vij = &v_cur[i * p * p];

				memcpy(Bij, xl_b, wXL * sizeof(double));
				/* 5) B_b := XR' * y */
				/*printf("DGEMV( T, %d, %d, %2f, %p, %d, %p, %d, %2f, %p, %d);\n",*/
				/*n, wXR, ONE, &x_cur[i * wXR * n], n, y_cur, iONE, ZERO, b_cur[i*p + wXL], iONE);*/
				dgemv_("T", 
						&n, &wXR, 
						&ONE, &x_cur[i * wXR * n], &n, y_cur, &iONE, 
						&ZERO, &Bij[wXL], &iONE);

				/* 5) W = XtZWinv * K^T */
				for( k = 0; k < wXL; k++ )
					dcopy_(&wXL, &xltxl[k*wXL], &iONE, &Vij[k*p], &iONE); // TL
				dgemm_("T", "N",
						&wXR, &wXL, &n, // m, n, k
						&ONE, &x_cur[i * wXR * n], &n, xl, &n, // KR KL'
						&ZERO, &Vij[wXL], &p); // xtSx BL
				dsyrk_("L", "T", 
						&wXR, &n, 
						&ONE, &x_cur[i * wXR * n], &n, 
						&ZERO, &Vij[wXL * p + wXL], &p);

				/* 8) W^-1 * y */
				dpotrf_(LOWER, &p, Vij, &p, &info);
				if (info != 0)
				{
					char err[STR_BUFFER_SIZE];
					snprintf(err, STR_BUFFER_SIZE, "dpotrf(V) failed (info: %d)", info);
					error_msg(err, 1);
				}
				dtrsv_(LOWER, NO_TRANS, NON_UNIT, &p, Vij, &p, Bij, &iONE);
				dtrsv_(LOWER,    TRANS, NON_UNIT, &p, Vij, &p, Bij, &iONE);

				// inv( X' inv(M) X)
				dpotri_(LOWER, &p, Vij, &p, &info);
				if (info != 0)
				{
					char err[STR_BUFFER_SIZE];
					snprintf(err, STR_BUFFER_SIZE, "dpotri failed (info: %d)", info);
					error_msg(err, 1);
				}
				//int p2 = p*p;
				//dscal_(&p2, &scorevar, xTSx, &iONE);
				for ( k = 0; k < p; k++ )
					Vij[k*p+k] = sqrt(Vij[k*p+k]);
			}

			if ( iter > 0)
			{
				fgls_aio_suspend( aiocb_b_prev_l, 1, NULL );
				fgls_aio_suspend( aiocb_v_prev_l, 1, NULL );
			}
			/* Write current B and V */
			struct aiocb *aiocb_b_cur_p;
			aiocb_b_cur_p = (struct aiocb *) aiocb_b_cur_l[0];
			fgls_aio_write( aiocb_b_cur_p,
							fileno( B_fp ), b_cur,
							x_inc * p * sizeof(double),
							(j * m * p + ib * p) * sizeof(double) );

			struct aiocb *aiocb_v_cur_p;
			aiocb_v_cur_p = (struct aiocb *) aiocb_v_cur_l[0];
			fgls_aio_write( aiocb_v_cur_p,
							fileno( V_fp ), v_cur,
							x_inc * p * p * sizeof(double),
							(j * m * p * p + ib * p * p) * sizeof(double) );

			swap_aiocb( &aiocb_x_cur_l, &aiocb_x_next_l );
			swap_buffers( &x_cur, &x_next);
			swap_aiocb( &aiocb_b_cur_l, &aiocb_b_prev_l );
			swap_buffers( &b_cur, &b_prev);
			swap_aiocb( &aiocb_v_cur_l, &aiocb_v_prev_l );
			swap_buffers( &v_cur, &v_prev);
			iter++;
		}
		/*gettimeofday(&t1, NULL);*/
		/*if (id == 0)*/
		/*{*/
		/*printf("Iter (%d - %d) time: %ld ms\n", jb, jb+y_inc, get_diff_ms(&t0, &t1));*/
		/*fflush(stdout);*/
		/*}*/
		swap_aiocb( &aiocb_y_cur_l, &aiocb_y_next_l );
		swap_buffers( &y_cur, &y_next);
	}

	fgls_aio_suspend( aiocb_x_cur_l, 1, NULL );
	fgls_aio_suspend( aiocb_y_cur_l, 1, NULL );
	fgls_aio_suspend( aiocb_b_prev_l, 1, NULL );
	fgls_aio_suspend( aiocb_v_prev_l, 1, NULL );

	/*pthread_exit(NULL);*/
	/*}*/

	fclose( XR_fp );
	fclose( Y_fp );
	fclose( B_fp );
	fclose( V_fp );

	for (i = 0; i < NUM_BUFFERS_PER_THREAD; i++) 
	{
		free( X[i] );
		free( Y[i] );
		free( B[i] );
		free( V[i] );
	}
	free( XL[0] );
	free( XL[1] );
	free( XL_b	);
	free( XLtXL );

	free( Phi );
	free( M );
	free( h );
	free( sigma );

  free( aiocb_x_cur_l  );
  free( aiocb_x_next_l );
  free( aiocb_y_cur_l  );
  free( aiocb_y_next_l );
  free( aiocb_b_prev_l );
  free( aiocb_b_cur_l  );
  free( aiocb_v_prev_l );
  free( aiocb_v_cur_l  );

	return 0;
}

