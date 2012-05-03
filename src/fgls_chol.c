#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <sys/types.h>
#include <sys/stat.h>

#include <sys/time.h>
#include <time.h>

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

/*
 * Cholesky-based solution of the Feasible Generalized Least-Squares problem
 */
int fgls_chol(int n, int p, int m, int t, int wXL, int wXR,
				 int x_b, int y_b, int num_threads, // y_b not used
				 char *Phi_path, char *h2_path, char *sigma2_path,
				 char *XL_path, char *XR_path, char *Y_path,
				 char *B_path, char *V_path)
{
	/* Configuration. Unnecessary (comes from the eigen-based implementation ) */
	FGLS_config_t cf;

	/* In-core operands */
	double *Phi;
	double *M;
	double *h;
	double *sigma;
	double alpha;
	double beta;

	/* sigma2.score */
	double *scoreB_t;
	double *scoreV_tl;
	double *scoreYmXB;
	double res_sigma;

	/* Out-of-core operands */
	double *X[NUM_BUFFERS_PER_THREAD];
	double *Y[NUM_BUFFERS_PER_THREAD];
	double *B[NUM_BUFFERS_PER_THREAD];
	double *V[NUM_BUFFERS_PER_THREAD];
	double *Bij, *Vij;

	/* Reusable data thanks to constant XL */
	double *XL;
	double *XL_orig; // XL and a copy (XL is overwritten at every iteration of j)
	double *B_t;  // Top part of b ( in inv(S) b )
	double *V_tl; // Top-Left part of V

	/* Data files */
	FILE *Phi_fp, *h_fp,  *sigma_fp,
		 *XL_fp,  *XR_fp, *Y_fp,
		 *B_fp,   *V_fp;


	/* BLAS / LAPACK constants */
	double ZERO = 0.0;
	double ONE = 1.0;
	double MINUS_ONE = -1.0;
	int iONE = 1;
	/* LAPACK error value */
	int info;

	/* iterators and auxiliar vars */
	int ib, i, j, k; // size_t
	int nn = n * n; // size_t
	char numths_str[STR_BUFFER_SIZE];

#if DEBUG
	/* Checking the input arguments */
	printf("n: %d\np: %d\nm: %d\nt: %d\nwXL: %d\nwXR: %d\n", n, p, m, t, wXL, wXR);
	printf("x_b: %d\ny_b: %d\nnths: %d\n", x_b, y_b, num_threads);
	printf("Phi: %s\nh2: %s\ns2: %s\n", Phi_path, h2_path, sigma2_path);
	printf("XL: %s\nXR: %s\nY: %s\n", XL_path, XR_path, Y_path);
	printf("B: %s\nV: %s\n", B_path, V_path);
#endif

	/* Fill in the config structure */
	initialize_config(
			&cf,
			n, p, m, t, wXL, wXR,
			x_b, 1, num_threads,
			Phi_path, h2_path, sigma2_path,
			XL_path, XR_path, Y_path,
			B_path, V_path
	);
	if ( y_b != 1 )
		fprintf(stderr, "[Warning] y_b not used (set to 1)\n");

	/* Memory allocation */
	// In-core
	Phi = ( double * ) fgls_malloc ( (size_t)cf.n * cf.n * sizeof(double) );
	M   = ( double * ) fgls_malloc ( (size_t)cf.n * cf.n * sizeof(double) );
	h     = ( double * ) fgls_malloc ( cf.t * sizeof(double) );
	sigma = ( double * ) fgls_malloc ( cf.t * sizeof(double) );

	XL_orig = ( double * ) fgls_malloc ( cf.wXL * cf.n * sizeof(double) );
	XL      = ( double * ) fgls_malloc ( cf.wXL * cf.n * sizeof(double) );
	B_t  = ( double * ) fgls_malloc ( cf.wXL * sizeof(double) );
	V_tl = ( double * ) fgls_malloc ( cf.wXL * cf.wXL * sizeof(double) );

	// sigma2.score
	scoreB_t  = ( double * ) fgls_malloc ( cf.wXL * sizeof(double) );
	scoreV_tl = ( double * ) fgls_malloc ( cf.wXL * cf.wXL * sizeof(double) );
	scoreYmXB = ( double * ) fgls_malloc ( cf.n * sizeof(double) );

	// Out-of-core 
	for (i = 0; i < NUM_BUFFERS_PER_THREAD; i++) 
	{
		X[i] = ( double * ) fgls_malloc ( (size_t) cf.x_b * cf.wXR * cf.n * sizeof(double) );
		/*X[i] = ( double * ) fgls_malloc ( (size_t) cf.x_b * cf.p * cf.n * sizeof(double) );*/
		Y[i] = ( double * ) fgls_malloc ( (size_t) cf.n * sizeof(double) );
		B[i] = ( double * ) fgls_malloc ( (size_t) cf.x_b * cf.p * sizeof(double) );
		V[i] = ( double * ) fgls_malloc ( (size_t) cf.x_b * cf.p * cf.p * sizeof(double) );
		/*V[i] = ( double * ) calloc ( cf.x_b * cf.p * cf.p, sizeof(double) );*/
	}

	/* Load in-core data */
	// Phi
	Phi_fp = fopen( cf.Phi_path, "rb" );
	sync_read( Phi, Phi_fp, cf.n * cf.n, 0 );
	fclose( Phi_fp );
	// h
	h_fp = fopen( cf.h_path, "r");
	sync_read(h, h_fp, cf.t, 0);
	fclose( h_fp );
	// sigma
	sigma_fp = fopen( cf.sigma_path, "r");
	sync_read(sigma, sigma_fp, cf.t, 0);
	fclose( sigma_fp );
	// XL: first column all ones, then columns from XL
	/*for ( i= 0; i < n; i++ )*/
	/*XL_orig[i] = 1.0;*/
	XL_fp = fopen( cf.XL_path, "rb" );
	sync_read( XL_orig, XL_fp, cf.wXL * cf.n, 0 );
	fclose( XL_fp );
	
	/* Files and pointers for out-of-core */
	XR_fp = fopen( cf.XR_path, "rb");
	Y_fp = fopen( cf.Y_path, "rb");
	B_fp = fopen( cf.B_path, "wb");
	V_fp = fopen( cf.V_path, "wb");

	double *x_cur  = X[0];
	double *x_next = X[1];
	double *y_cur  = Y[0];
	double *y_next = Y[1];
	double *b_cur  = B[0];
	double *b_prev = B[1];
	double *v_cur  = V[0];
	double *v_prev = V[1];

	/* Asynchronous IO data structures */
	struct aiocb aiocb_x_cur,  aiocb_x_next,
				 aiocb_y_cur,  aiocb_y_next,
				 aiocb_b_prev, aiocb_b_cur,
				 aiocb_v_prev, aiocb_v_cur;

	const struct aiocb ** aiocb_x_cur_l, // = { &aiocb_x_cur }, 
		               ** aiocb_x_next_l,// = { &aiocb_x_next },
				       ** aiocb_y_cur_l, // = { &aiocb_y_cur },
				       ** aiocb_y_next_l,// = { &aiocb_y_next },
				       ** aiocb_b_prev_l,// = { aiocb_b_prev },
				       ** aiocb_b_cur_l, // = { aiocb_b_cur };
				       ** aiocb_v_prev_l,// = { aiocb_v_prev },
				       ** aiocb_v_cur_l; // = { aiocb_v_cur };

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

#if VAMPIR
      VT_USER_START("READ_X");
#endif
	/* Read first block of XR's */
	fgls_aio_read( &aiocb_x_cur,
			       fileno( XR_fp ), x_cur,
				   (size_t)MIN( x_b, m ) * wXR * n * sizeof(double), 0 );
#if VAMPIR
      VT_USER_END("READ_X");
#endif
#if VAMPIR
      VT_USER_START("READ_Y");
#endif
	/* Read first Y */
	fgls_aio_read( &aiocb_y_cur,
			       fileno( Y_fp ), y_cur,
				   (size_t)n * sizeof(double), 0 );
#if VAMPIR
      VT_USER_END("READ_Y");
#endif

	int iter = 0;
	for ( j = 0; j < t; j++ )
	{
		/* Set the number of threads for the multi-threaded BLAS */
		snprintf(numths_str, STR_BUFFER_SIZE, "%d", cf.NUM_COMPUTE_THREADS);
#if defined GOTO
		setenv("GOTO_NUM_THREADS", numths_str, 1);
#elif defined MKL
		setenv("MKL_NUM_THREADS", numths_str, 1);
#else
		setenv("OMP_NUM_THREADS", numths_str, 1);
#endif

#if VAMPIR
      VT_USER_START("READ_Y");
#endif
		/* Read next Y */
		struct aiocb *aiocb_y_next_p = (struct aiocb *)aiocb_y_next_l[0];
		fgls_aio_read( aiocb_y_next_p,
					   fileno( Y_fp ), y_next,
					   j + 1 >= t ? 0 : (size_t)n * sizeof(double),
					   j + 1 >= t ? 0 : (off_t)(j+1) * n * sizeof(double) );
#if VAMPIR
      VT_USER_END("READ_Y");
#endif

#if VAMPIR
      VT_USER_START("COMP_J");
#endif
		/* M := sigma * ( h^2 Phi - (1 - h^2) I ) */
		memcpy( M, Phi, (size_t)n * n * sizeof(double) );
		alpha = h[j] * sigma[j];
		beta  = (1 - h[j]) * sigma[j];
		dscal_(&nn, &alpha, M, &iONE);
		for ( i = 0; i < n; i++ )
			M[i*n + i] = M[i*n + i] + beta;


		/* L * L' = M */
		/*struct timeval t0, t1;*/
		/*gettimeofday(&t0, NULL);*/
		dpotrf_(LOWER, &n, M, &n, &info);
		/*gettimeofday(&t1, NULL);*/
		/*printf("Chol:  %ld msecs\n", get_diff_ms( &t0, &t1 ));*/
		if (info != 0)
		{
			char err[STR_BUFFER_SIZE];
			snprintf(err, STR_BUFFER_SIZE, "dpotrf(M) failed (info: %d)", info);
			error_msg(err, 1);
		}

		/* XL := inv(L) XL */
		memcpy( XL, XL_orig, wXL * n * sizeof(double) );
		dtrsm_(LEFT, LOWER, NO_TRANS, NON_UNIT, &n, &wXL, &ONE, M, &n, XL, &n);

#if VAMPIR
      VT_USER_START("WAIT_Y");
#endif
		/* Wait until current Y is available */
	  /*printf("Waiting for %lu bytes read (Y)\n", aiocb_y_cur.aio_nbytes);*/
		fgls_aio_suspend( aiocb_y_cur_l, 1, NULL );
#if VAMPIR
      VT_USER_END("WAIT_Y");
#endif
	    // Copy y for sigma2.score
		for( k = 0; k < n; k++ )
			scoreYmXB[k] = y_cur[k];

		/* y := inv(L) y */
		dtrsv_(LOWER, NO_TRANS, NON_UNIT, &n, M, &n, y_cur, &iONE);

		/* B_t := XL' * y */
		dgemv_(TRANS, &n, &wXL, &ONE, XL, &n, y_cur, &iONE, &ZERO, B_t, &iONE);

		/* V_tl := XL' * XL */
		dsyrk_(LOWER, TRANS, &wXL, &n, &ONE, XL, &n, &ZERO, V_tl, &wXL);
#if VAMPIR
      VT_USER_END("COMP_J");
#endif
		/* Compute sigma2.score */
		// copy B_t and V_tl
		for( k = 0; k < wXL; k++ )
			scoreB_t[k] = B_t[k];
		for( k = 0; k < wXL*wXL; k++ )
			scoreV_tl[k] = V_tl[k];
		// XB
		dpotrf_(LOWER, &wXL, scoreV_tl, &wXL, &info);
		if (info != 0)
		{
			char err[STR_BUFFER_SIZE];
			snprintf(err, STR_BUFFER_SIZE, "dpotrf(scoreV) failed (info: %d) - i: %d", info, i);
			error_msg(err, 1);
		}
		dtrsv_(LOWER, NO_TRANS, NON_UNIT, &wXL, scoreV_tl, &wXL, scoreB_t, &iONE);
		dtrsv_(LOWER,    TRANS, NON_UNIT, &wXL, scoreV_tl, &wXL, scoreB_t, &iONE);
		// YmXB
		dgemv_(NO_TRANS, 
				&n, &wXL, 
				&MINUS_ONE, XL_orig, &n, scoreB_t, &iONE,
				&ONE, scoreYmXB, &iONE);
		// residual sigma
		dtrsv_(LOWER, NO_TRANS, NON_UNIT, 
				&n, M, &n, scoreYmXB, &iONE);
		res_sigma = ddot_(&n, scoreYmXB, &iONE, scoreYmXB, &iONE) / (n - wXL);

		for (ib = 0; ib < m; ib += x_b) 
		{
#if VAMPIR
      VT_USER_START("READ_X");
#endif
			/* Read next block of XR's */
			struct aiocb *aiocb_x_next_p = (struct aiocb *)aiocb_x_next_l[0];
			fgls_aio_read( aiocb_x_next_p,
					       fileno( XR_fp ), x_next,
						   ((size_t)ib + x_b) >= m ? (size_t)MIN( x_b, m ) * wXR * n * sizeof(double) : MIN((size_t)x_b, m - ((size_t)ib + x_b)) * wXR * n * sizeof(double),
						   ((off_t)ib + x_b) >= m ? 0 : ((off_t)ib + x_b) * wXR * n * sizeof(double) );
#if VAMPIR
      VT_USER_END("READ_X");
#endif

#if VAMPIR
      VT_USER_START("WAIT_X");
#endif
			/* Wait until current block of XR's is available */
	  /*printf("Waiting for %lu bytes read (X)\n", aiocb_x_cur.aio_nbytes);*/
			fgls_aio_suspend( aiocb_x_cur_l, 1, NULL );
#if VAMPIR
      VT_USER_END("WAIT_X");
#endif

		/* Set the number of threads for the multi-threaded BLAS */
		snprintf(numths_str, STR_BUFFER_SIZE, "%d", cf.NUM_COMPUTE_THREADS);
#if defined GOTO
		setenv("GOTO_NUM_THREADS", numths_str, 1);
#elif defined MKL
		setenv("MKL_NUM_THREADS", numths_str, 1);
#else
		setenv("OMP_NUM_THREADS", numths_str, 1);
#endif

#if VAMPIR
      VT_USER_START("COMP_IB");
#endif
			/* XR := inv(L) XR */
			int x_inc = MIN(x_b, m - ib);
			int rhss  = wXR * x_inc;
			dtrsm_(LEFT, LOWER, NO_TRANS, NON_UNIT, &n, &rhss, &ONE, M, &n, x_cur, &n);
#if VAMPIR
      VT_USER_END("COMP_IB");
#endif

			/* Set the number of threads for the multi-threaded BLAS to 1.
			 * The innermost loop is parallelized using OPENMP */
		#if defined GOTO
			setenv("GOTO_NUM_THREADS", "1", 1);
		#elif defined MKL
			setenv("MKL_NUM_THREADS",  "1", 1);
		#else
			setenv("OMP_NUM_THREADS",  "1", 1);
		#endif
#if VAMPIR
      VT_USER_START("COMP_I");
#endif
			#pragma omp parallel for private(Bij, Vij, i, k, info) schedule(static) num_threads(cf.NUM_COMPUTE_THREADS)
			for (i = 0; i < x_inc; i++)
			{
				Bij = &b_cur[i * p];
				Vij = &v_cur[i * p * p];

				/* Building B */
				// Copy B_T
				memcpy(Bij, B_t, wXL * sizeof(double));
				// B_B := XR' * y
				dgemv_("T", 
						&n, &wXR, 
						&ONE, &x_cur[i * wXR * n], &n, y_cur, &iONE, 
						&ZERO, &Bij[wXL], &iONE);

				/* Building V */
				// Copy V_TL
				for( k = 0; k < wXL; k++ )
					dcopy_(&wXL, &V_tl[k*wXL], &iONE, &Vij[k*p], &iONE); // V_TL
				// V_BL := XR' * XL
				dgemm_("T", "N",
						&wXR, &wXL, &n,
						&ONE, &x_cur[i * wXR * n], &n, XL, &n,
						&ZERO, &Vij[wXL], &p); // V_BL
				// V_BR := XR' * XR
				dsyrk_("L", "T", 
						&wXR, &n, 
						&ONE, &x_cur[i * wXR * n], &n, 
						&ZERO, &Vij[wXL * p + wXL], &p); // V_BR

				/* B := inv(V) * B */
				dpotrf_(LOWER, &p, Vij, &p, &info);
				if (info != 0)
				{
					char err[STR_BUFFER_SIZE];
					snprintf(err, STR_BUFFER_SIZE, "dpotrf(V) failed (info: %d) - i: %d", info, i);
					error_msg(err, 1);
				}
				dtrsv_(LOWER, NO_TRANS, NON_UNIT, &p, Vij, &p, Bij, &iONE);
				dtrsv_(LOWER,    TRANS, NON_UNIT, &p, Vij, &p, Bij, &iONE);

				/* V := inv( X' inv(M) X) */
				dpotri_(LOWER, &p, Vij, &p, &info);
				if (info != 0)
				{
					char err[STR_BUFFER_SIZE];
					snprintf(err, STR_BUFFER_SIZE, "dpotri failed (info: %d)", info);
					error_msg(err, 1);
				}
				int p2 = p*p;
				dscal_(&p2, &res_sigma, Vij, &iONE);
				for ( k = 0; k < p; k++ )
					Vij[k*p+k] = sqrt(Vij[k*p+k]);
			}
#if VAMPIR
      VT_USER_END("COMP_I");
#endif

#if VAMPIR
      VT_USER_START("WAIT_BV");
#endif
			/* Wait until the previous blocks of B's and V's are written */
			if ( iter > 0)
			{
				/*printf("Waiting for %lu bytes written (b)\n", aiocb_b_prev.aio_nbytes);*/
				fgls_aio_suspend( aiocb_b_prev_l, 1, NULL );
				/*printf("Waiting for %lu bytes written (v)\n", aiocb_v_prev.aio_nbytes);*/
				fgls_aio_suspend( aiocb_v_prev_l, 1, NULL );
			}
#if VAMPIR
      VT_USER_END("WAIT_BV");
#endif
			/* Write current blocks of B's and V's */
#if VAMPIR
      VT_USER_START("WRITE_BV");
#endif
			struct aiocb *aiocb_b_cur_p;
			aiocb_b_cur_p = (struct aiocb *) aiocb_b_cur_l[0];
			fgls_aio_write( aiocb_b_cur_p,
							fileno( B_fp ), b_cur,
							(size_t)x_inc * p * sizeof(double),
							((off_t)j * m * p + (off_t)ib * p) * sizeof(double) );

			struct aiocb *aiocb_v_cur_p;
			aiocb_v_cur_p = (struct aiocb *) aiocb_v_cur_l[0];
			fgls_aio_write( aiocb_v_cur_p,
							fileno( V_fp ), v_cur,
							(size_t)x_inc * p * p * sizeof(double),
							((off_t)j * m * p * p + (off_t)ib * p * p) * sizeof(double) );
#if VAMPIR
      VT_USER_END("WRITE_BV");
#endif

			/* Swap buffers */
			swap_aiocb( &aiocb_x_cur_l, &aiocb_x_next_l );
			swap_buffers( &x_cur, &x_next);
			swap_aiocb( &aiocb_b_cur_l, &aiocb_b_prev_l );
			swap_buffers( &b_cur, &b_prev);
			swap_aiocb( &aiocb_v_cur_l, &aiocb_v_prev_l );
			swap_buffers( &v_cur, &v_prev);
			iter++;
		}
		/* Swap buffers */
		swap_aiocb( &aiocb_y_cur_l, &aiocb_y_next_l );
		swap_buffers( &y_cur, &y_next);
	}

#if VAMPIR
      VT_USER_START("WAIT_ALL");
#endif
	/* Wait for the remaining IO operations issued */
	  /*printf("Last Waiting");*/
	fgls_aio_suspend( aiocb_x_cur_l, 1, NULL );
	fgls_aio_suspend( aiocb_y_cur_l, 1, NULL );
	fgls_aio_suspend( aiocb_b_prev_l, 1, NULL );
	fgls_aio_suspend( aiocb_v_prev_l, 1, NULL );
#if VAMPIR
      VT_USER_END("WAIT_ALL");
#endif

	/* Clean-up */
	fclose( XR_fp );
	fclose( Y_fp );
	fclose( B_fp );
	fclose( V_fp );

	free( Phi );
	free( M );
	free( h );
	free( sigma );

	free( XL );
	free( XL_orig );
	free( B_t  );
	free( V_tl );

	free( scoreB_t );
	free( scoreV_tl );
	free( scoreYmXB );

	for (i = 0; i < NUM_BUFFERS_PER_THREAD; i++) 
	{
		free( X[i] );
		free( Y[i] );
		free( B[i] );
		free( V[i] );
	}

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

