#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <sys/types.h>
#include <sys/stat.h>

#include <sys/time.h>
#include <time.h>

#include <pthread.h>
#include <aio.h>

#include "blas.h"
#include "lapack.h"
#include "options.h"
#include "common.h"
#include "io.h"
#include "timing.h"
#include "fgls_eigen.h"

#if VAMPIR
  #include "vt_user.h"
#endif

#include <stdint.h>

int    preloop(FGLS_config_t *cf, double *Phi, double *Z, double *W, double *res_sigma);
void   eigenDec(int n, double *Phi, double *Z, double *W);
void * ooc_gemm( void *in );
void residual_sigma( FGLS_config_t *fgls_cf, ooc_res_sigma_t *cf ); 
void * ooc_loops(void* in);

/*
 * Eigen-based computation of the Feasible Generalized Least-Squares problem
 */
int fgls_eigen(int n, int p, int m, int t, int wXL, int wXR,
               int x_b, int y_b, int num_threads,
               char *Phi_path, char *h2_path, char *sigma2_path,
               char *XL_path, char *XR_path, char *Y_path,
               char *B_path, char *V_path)
{
	/* Problem configuration */
	FGLS_config_t cf;

	/* Preloop computation */
	double *Phi,
		   *Z;
	FILE *Phi_fp, *h_fp, *sigma_fp,
	     *XL_fp;

	/* Loops computation */
	// Thread handling 
	int rc;
	void *retval;
	pthread_t *compute_threads;
	// Data structures
	ooc_loops_t loops_t;
	ooc_loops_t *loops_t_comp;

	// iterators and auxiliar vars
	int i;
	char numths_str[STR_BUFFER_SIZE];

	/* Check input values */
	/*printf("n: %d\np: %d\nm: %d\nt: %d\nwXL: %d\nwXR: %d\n", n, p, m, t, wXL, wXR);*/
	/*printf("x_b: %d\ny_b: %d\nnths: %d\n", x_b, y_b, num_threads);*/
	/*printf("Phi: %s\nh2: %s\ns2: %s\n", cf.Phi_path, cf.h_path, cf.sigma_path);*/
	/*printf("XL: %s\nXR: %s\nY: %s\n", cf.XL_path, cf.XR_path, cf.Y_path);*/
	/*printf("B: %s\nV: %s\n", cf.B_path, cf.V_path);*/

	/* Fill in the config structure */
	initialize_config(
		&cf,
		n, p, m, t, wXL, wXR,
		x_b, y_b, num_threads,
		Phi_path, h2_path, sigma2_path,
		XL_path, XR_path, Y_path,
		B_path, V_path
	);


#if VAMPIR
	VT_USER_START("PRELOOP");
#endif
	/* Allocate memory for the eigendecomposition */
	Phi = ( double * ) fgls_malloc ( cf.n * cf.n * sizeof(double) );
	Z = ( double * ) fgls_malloc ( cf.n * cf.n * sizeof(double) );
	loops_t.W = ( double * ) fgls_malloc ( cf.n * sizeof(double) );
	loops_t.res_sigma = ( double * ) fgls_malloc ( cf.t * sizeof(double) );
	
	/* Load Phi */
	Phi_fp = fopen( cf.Phi_path, "rb" );
	sync_read( Phi, Phi_fp, cf.n * cf.n, 0 );
	fclose( Phi_fp );
	
	/* Compute the pre-loop operations */
	snprintf(numths_str, STR_BUFFER_SIZE, "%d", cf.NUM_COMPUTE_THREADS);
#if defined GOTO
	setenv("GOTO_NUM_THREADS", numths_str, 1);
#elif defined MKL
	setenv("MKL_NUM_THREADS", numths_str, 1);
#else
	setenv("OMP_NUM_THREADS", numths_str, 1);
#endif
	preloop(&cf, Phi, Z, loops_t.W, loops_t.res_sigma);
#if defined GOTO
	setenv("GOTO_NUM_THREADS", "1", 1);
#elif defined MKL
	setenv("MKL_NUM_THREADS", "1", 1);
#else
	setenv("OMP_NUM_THREADS", "1", 1);
#endif

	/* Clean-up */
	free( Phi );
	free( Z );
#if VAMPIR
	VT_USER_END("PRELOOP");
#endif

#if VAMPIR
	VT_USER_START("LOOPS");
#endif
	loops_t.cf = &cf;
	/* Allocating memory for the computation of the double loop */
	// In-core data
	loops_t.h     = ( double * ) fgls_malloc ( cf.t * sizeof(double) );
	loops_t.sigma = ( double * ) fgls_malloc ( cf.t * sizeof(double) );
	loops_t.alpha = ( double * ) fgls_malloc ( cf.NUM_COMPUTE_THREADS * cf.y_b * sizeof(double) );
	loops_t.beta  = ( double * ) fgls_malloc ( cf.NUM_COMPUTE_THREADS * cf.y_b * sizeof(double) );
	loops_t.Winv  = ( double * ) fgls_malloc ( cf.NUM_COMPUTE_THREADS * cf.x_b * cf.n * sizeof(double) );
	// The original in XL[0]
	// The copy to overwrite in XL[1]
	loops_t.XL_orig = ( double * ) fgls_malloc ( cf.wXL * cf.n * sizeof(double) );
	loops_t.XL[0] = ( double * ) fgls_malloc ( cf.wXL * cf.n * sizeof(double) );
	loops_t.XL[1] = ( double * ) fgls_malloc ( cf.NUM_COMPUTE_THREADS * cf.wXL * cf.n * cf.y_b * sizeof(double) );
	// Reusable B_Top and V_TopLeft
	loops_t.B_t  = ( double * ) fgls_malloc ( cf.NUM_COMPUTE_THREADS * cf.wXL * cf.y_b * sizeof(double) );
	loops_t.V_tl = ( double * ) fgls_malloc ( cf.NUM_COMPUTE_THREADS * cf.wXL * cf.wXL * cf.y_b * sizeof(double) );

	// Double-buffering for the operands that do not fit in RAM
	for (i = 0; i < NUM_BUFFERS_PER_THREAD; i++) 
	{
		loops_t.X[i] = ( double* ) fgls_malloc ( cf.NUM_COMPUTE_THREADS * cf.x_b * cf.p * cf.n * sizeof(double) );
		loops_t.Y[i] = ( double* ) fgls_malloc ( cf.NUM_COMPUTE_THREADS * cf.y_b * cf.n * sizeof(double) );
		loops_t.B[i] = ( double* ) fgls_malloc ( cf.NUM_COMPUTE_THREADS * cf.x_b * cf.y_b * cf.p * sizeof(double) );
		loops_t.V[i] = ( double* ) fgls_malloc ( cf.NUM_COMPUTE_THREADS * cf.x_b * cf.y_b * cf.p * cf.p * sizeof(double) );
	}

	/* Load in-core data */
	// h
	h_fp = fopen( cf.h_path, "r");
	sync_read(loops_t.h, h_fp, cf.t, 0);
	fclose( h_fp );
	// sigma
	sigma_fp = fopen( cf.sigma_path, "r");
	sync_read(loops_t.sigma, sigma_fp, cf.t, 0);
	fclose( sigma_fp );
	// XL: first column all ones, then columns from XL
	XL_fp = fopen( cf.ZtXL_path, "rb" );
	sync_read( loops_t.XL[0], XL_fp, cf.wXL * cf.n, 0 );
	fclose( XL_fp );
	XL_fp = fopen( cf.XL_path, "rb" );
	sync_read( loops_t.XL_orig, XL_fp, cf.wXL * cf.n, 0 );
	fclose( XL_fp );

	/* Files for out-of-core data: XR, Y, B and V */
	loops_t.XR_fp = fopen( cf.ZtXR_path, "rb");
	loops_t.Y_fp  = fopen( cf.ZtY_path, "rb");
	loops_t.B_fp  = fopen( cf.B_path, "wb");
	loops_t.V_fp  = fopen( cf.V_path, "wb");

	/*printf("LOOPS\n");*/

	/* Spawn threads */
	loops_t_comp = ( ooc_loops_t* ) fgls_malloc ( cf.NUM_COMPUTE_THREADS * sizeof(ooc_loops_t) );
	compute_threads = ( pthread_t * ) fgls_malloc ( cf.NUM_COMPUTE_THREADS * sizeof(pthread_t) );
	for (i = 0; i < cf.NUM_COMPUTE_THREADS; i++) 
	{
		memcpy((void*)&loops_t_comp[i], (void*)&loops_t, sizeof(ooc_loops_t));
		loops_t_comp[i].id = i;

		rc = pthread_create( &compute_threads[i], NULL, ooc_loops, (void*)&loops_t_comp[i] );
		if (rc) 
		{
			char err[STR_BUFFER_SIZE];
			snprintf(err, STR_BUFFER_SIZE, "pthread_create() returned %d\n", rc);
			error_msg(err, 1);
		}
	}

	/* Wait for the spawned threads */
	for (i = 0; i < cf.NUM_COMPUTE_THREADS; i++)
		pthread_join(compute_threads[i], &retval);
#if VAMPIR
	VT_USER_START("END");
#endif

	/*printf("LOOPS Done\n");*/

	/* Clean-up */
	fclose( loops_t.XR_fp );
	fclose( loops_t.Y_fp );
	fclose( loops_t.B_fp );
	fclose( loops_t.V_fp );

	for (i = 0; i < NUM_BUFFERS_PER_THREAD; i++) 
	{
		free( loops_t.X[i] );
		free( loops_t.Y[i] );
		free( loops_t.B[i] );
		free( loops_t.V[i] );
	}
	free( loops_t.XL_orig );
	free( loops_t.XL[0] );
	free( loops_t.XL[1] );
	free( loops_t.B_t  );
	free( loops_t.V_tl );

	free( loops_t.h );
	free( loops_t.sigma );
	free( loops_t.res_sigma );
	free( loops_t.W );

	free( loops_t.alpha );
	free( loops_t.beta );
	free( loops_t.Winv );

	free( loops_t_comp );
	free( compute_threads );

	return 0;
}

void* ooc_loops(void* in) 
{
  ooc_loops_t *loops_t = ( ooc_loops_t* ) in;
  FGLS_config_t *cf = loops_t->cf;

  /* Dimensions of the problem */
  int m = cf->m,
	  n = cf->n,
	  p = cf->p,
	  t = cf->t,
	  wXL = cf->wXL,
	  wXR = cf->wXR,
	  x_b = cf->x_b,
	  y_b = cf->y_b,
	  id = loops_t->id;

  /* In-core data */
  double *W      =  loops_t->W;
  double *h      =  loops_t->h;
  double *sigma  =  loops_t->sigma;
  double *res_sigma  =  loops_t->res_sigma;
  double *alpha  = loops_t->alpha, 
		 *beta   = loops_t->beta;
  double *Winv = &loops_t->Winv[id * y_b * n];

  /* Reusable data (common XL) */
  double *XL     = &loops_t->XL[1][id * wXL * n * y_b];
  double *B_t    = &loops_t->B_t[id * wXL * y_b];
  double *V_tl   = &loops_t->V_tl[id * wXL * wXL * y_b];

	/* sigma2.score */
  /*double *scoreB_t;*/
  /*double *scoreV_tl;*/
  /*double *scoreYmXB;*/
  /*double *res_sigma;*/

  /* Double buffering pointers */
  double *x_cur  = &loops_t->X[0][id * n * wXR * x_b];
  double *x_next = &loops_t->X[1][id * n * wXR * x_b];
  double *y_cur  = &loops_t->Y[0][id * y_b * n];
  double *y_next = &loops_t->Y[1][id * y_b * n];
  double *b_cur  = &loops_t->B[0][id * p * x_b * y_b];
  double *b_prev = &loops_t->B[1][id * p * x_b * y_b];
  double *v_cur  = &loops_t->V[0][id * p * p * x_b * y_b];
  double *v_prev = &loops_t->V[1][id * p * p * x_b * y_b];

  /* BLAS / LAPACK constants */
  double ZERO = 0.0;
  double ONE = 1.0;
  /*double MINUS_ONE = -1.0;*/
  int iONE = 1;
  /* LAPACK error value */
  int info;

  /* iterators and auxiliar vars */
  int ib, jb, i, j, k, l, ll;
  int x_inc, y_inc;
  double *Bij, *Vij;

  /* XR used many times and overwritten, need to work on a copy */
  double *x_copy = (double *) fgls_malloc (x_b * wXR * n * sizeof(double));

  /* sigma2.score */
  /*scoreB_t  = (double *) fgls_malloc (wXL * y_b * sizeof(double));*/
  /*scoreV_tl = (double *) fgls_malloc (wXL * wXL * y_b * sizeof(double));*/
  /*scoreYmXB = (double *) fgls_malloc (n * y_b * sizeof(double));*/
  /*res_sigma = (double *) fgls_malloc (y_b * sizeof(double));*/

  /* Asynchronous IO data structures */
  struct aiocb  aiocb_x_cur,   aiocb_x_next,
			    aiocb_y_cur,   aiocb_y_next,
			   *aiocb_b_prev, *aiocb_b_cur,
			   *aiocb_v_prev, *aiocb_v_cur;
  aiocb_b_prev = (struct aiocb *) fgls_malloc (y_b * sizeof(struct aiocb));
  aiocb_b_cur  = (struct aiocb *) fgls_malloc (y_b * sizeof(struct aiocb));
  aiocb_v_prev = (struct aiocb *) fgls_malloc (y_b * sizeof(struct aiocb));
  aiocb_v_cur  = (struct aiocb *) fgls_malloc (y_b * sizeof(struct aiocb));

  struct aiocb const ** aiocb_x_cur_l,  //[1]  = { &aiocb_x_cur }, 
			         ** aiocb_x_next_l, //[1] = { &aiocb_x_next },
			         ** aiocb_y_cur_l,  //[1]  = { &aiocb_y_cur },
			         ** aiocb_y_next_l, //[1] = { &aiocb_y_next },
			         ** aiocb_b_prev_l, //[1] = {  aiocb_b_prev },
			         ** aiocb_b_cur_l,  //[1]  = {  aiocb_b_cur },
			         ** aiocb_v_prev_l, //[1] = {  aiocb_v_prev },
			         ** aiocb_v_cur_l;  //[1]  = {  aiocb_v_cur };

  aiocb_x_cur_l  = (const struct aiocb **) fgls_malloc (sizeof(struct aiocb *));
  aiocb_x_next_l = (const struct aiocb **) fgls_malloc (sizeof(struct aiocb *));
  aiocb_y_cur_l  = (const struct aiocb **) fgls_malloc (sizeof(struct aiocb *));
  aiocb_y_next_l = (const struct aiocb **) fgls_malloc (sizeof(struct aiocb *));
  aiocb_b_prev_l = (const struct aiocb **) fgls_malloc (y_b * sizeof(struct aiocb *));
  aiocb_b_cur_l  = (const struct aiocb **) fgls_malloc (y_b * sizeof(struct aiocb *));
  aiocb_v_prev_l = (const struct aiocb **) fgls_malloc (y_b * sizeof(struct aiocb *));
  aiocb_v_cur_l  = (const struct aiocb **) fgls_malloc (y_b * sizeof(struct aiocb *));

  aiocb_x_cur_l[0]  = &aiocb_x_cur;
  aiocb_x_next_l[0] = &aiocb_x_next;
  aiocb_y_cur_l[0]  = &aiocb_y_cur;
  aiocb_y_next_l[0] = &aiocb_y_next;
  for ( i = 0; i < y_b; i++ )
  {
	  aiocb_b_prev_l[i] = &aiocb_b_prev[i];
	  aiocb_b_cur_l[i]  = &aiocb_b_cur[i];
	  aiocb_v_prev_l[i] = &aiocb_v_prev[i];
	  aiocb_v_cur_l[i]  = &aiocb_v_cur[i];
  }

  /* Read initial XR's and Y */
  fgls_aio_read( &aiocb_x_cur, 
		         fileno( loops_t->XR_fp ), x_cur,
				 (size_t)MIN( x_b, m ) * wXR * n * sizeof(double), 0 );

  fgls_aio_read( &aiocb_y_cur,
		         fileno( loops_t->Y_fp ), y_cur,
				 id * y_b >= t ? 0 : (size_t)n * MIN(y_b, t - id * y_b) * sizeof(double),
				 id * y_b >= t ? 0 : (off_t)id * n * y_b * sizeof(double) );

  int iter = 0;
  struct timeval t0, t1;
  for (jb = id * y_b; jb < t; jb += cf->NUM_COMPUTE_THREADS * y_b) 
  {
	  /*printf("Iter jb: %d\n", jb);*/
    gettimeofday(&t0, NULL);
	y_inc = MIN( y_b, t - jb );
#if VAMPIR
      VT_USER_START("READ_Y");
#endif
	/* Read next Y */
	struct aiocb *aiocb_y_next_p = (struct aiocb *)aiocb_y_next_l[0];
	fgls_aio_read( aiocb_y_next_p,
			       fileno( loops_t->Y_fp ), y_next,
				   jb + cf->NUM_COMPUTE_THREADS * y_b >= t ? 0 : (size_t)n * MIN(y_b, t - (jb + cf->NUM_COMPUTE_THREADS * y_b)) * sizeof(double),
				   jb + cf->NUM_COMPUTE_THREADS * y_b >= t ? 0 : (off_t)(jb + y_b * cf->NUM_COMPUTE_THREADS) * n * sizeof(double) );

#if VAMPIR
      VT_USER_END("READ_Y");
#endif

#if VAMPIR
      VT_USER_START("WAIT_X");
#endif
	/* Copy XL */
	for (ll = 0; ll < y_inc; ll++)
		memcpy( &XL[ll * wXL * n], loops_t->XL[0], wXL * n * sizeof(double) );
#if VAMPIR
      VT_USER_END("WAIT_X");
#endif
#if VAMPIR
    VT_USER_START("WAIT_Y");
#endif
	/* Wait until the current Y is available */
	fgls_aio_suspend( aiocb_y_cur_l, 1, NULL );
#if VAMPIR
    VT_USER_END("WAIT_Y");
#endif
#if VAMPIR
    VT_USER_START("COMP_LOOP_Y_CODE");
#endif
	/* Copy y for sigma2.score */
	/*memcpy( scoreYmXB, y_cur, n * y_inc * sizeof(double) );*/

	/* Set the scalars alpha and beta to compute: 
	 *     Winv = alpha W + beta I */
    for (k = 0; k < y_inc; k++)
	{
    	alpha[k] = sigma[jb + k] * h[jb + k]; // * h[j];
    	beta[k]  = sigma[jb + k] * (1 - h[jb + k]); // * h[j]);
	}

    /* Winv := sqrt( inv( alpha W - beta I ) ) */
    // Best order? sqrt - inv
	//
	// Possibly GER
    for (k = 0; k < y_inc; k++)
        for (l = 0; l < n; l++)
            Winv[k*n + l] = sqrt(1.0 / (alpha[k] * W[l] + beta[k]));

    /* y := sqrt(Winv) * Z' * y */
    for (k = 0; k < y_inc; k++)
        for (l = 0; l < n; l++)
            y_cur[k*n+l] *= Winv[k*n+l];

    /* XL := sqrt(Winv) * Z' * XL */
	for (ll = 0; ll < y_inc; ll++)
	  for ( k = 0; k < wXL; k++ )
		  for ( l = 0; l < n; l++ )
			  XL[ ll * wXL * n + k * n + l ] *= Winv[ll * n + l];
      
    /* B_t := XL' * XL */
	for (ll = 0; ll < y_inc; ll++)
      dgemv_("T", &n, &wXL, &ONE, &XL[ll * wXL * n], &n, &y_cur[ll * n], &iONE, &ZERO, &B_t[ll * wXL], &iONE);

	/* V_tl := Compute Top-left part of V */
	for (ll = 0; ll < y_inc; ll++)
	{
		dsyrk_("L", "T", // LOWER, NO_TRANS, 
				&wXL, &n, // n, k
				&ONE, &XL[ll * wXL * n], &n, // KL KL' 
				&ZERO, &V_tl[ll * wXL * wXL], &wXL); // V_TL
	}

	/* Compute sigma2.score's */
	// copy B_t and V_tl
	/*memcpy( scoreB_t,  B_t,  y_inc * wXL * sizeof(double) );*/
	/*memcpy( scoreV_tl, V_tl, y_inc * wXL * wXL * sizeof(double) );*/
	// XB
	/*double *sc_B_t, *sc_V_tl;*/
	/*for ( ll = 0; ll < y_inc; ll++ )*/
	/*{*/
	/*sc_B_t  = &scoreB_t [ ll * wXL ];*/
	/*sc_V_tl = &scoreV_tl[ ll * wXL * wXL ];*/
	/**/
	/*dpotrf_(LOWER, &wXL, sc_V_tl, &wXL, &info);*/
	/*if (info != 0)*/
	/*{*/
	/*char err[STR_BUFFER_SIZE];*/
	/*snprintf(err, STR_BUFFER_SIZE, "sigma2.score: dpotrf(scoreV(%d)) failed (info: %d) - i: %d", ll, info, i);*/
	/*error_msg(err, 1);*/
	/*}*/
	/*printf("V_tl[%d] = %.16e\n", ll, sc_V_tl[0]);*/
	/*dtrsv_(LOWER, NO_TRANS, NON_UNIT, &wXL, sc_V_tl, &wXL, sc_B_t, &iONE);*/
	/*dtrsv_(LOWER,    TRANS, NON_UNIT, &wXL, sc_V_tl, &wXL, sc_B_t, &iONE);*/
	/*printf("B_t [%d] = %.16e\n", ll, sc_B_t[0]);*/
	/*}*/
	// YmXB
	/*double *XL_orig = XL[0];*/
	/*printf("XL[0] = %.16e\n", loops_t->XL_orig[0]);*/
	/*printf("YmXB[0] = %.16e\n", scoreYmXB[0]);*/
	/*printf("YmXB[1] = %.16e\n", scoreYmXB[n]);*/
	/*dgemm_(NO_TRANS, NO_TRANS,*/
	/*&n, &y_inc, &wXL, */
	/*&MINUS_ONE, loops_t->XL_orig, &n, scoreB_t, &wXL,*/
	/*&ONE, scoreYmXB, &n);*/
	/*printf("YmXB[0] = %.16e\n", scoreYmXB[0]);*/
	/*printf("YmXB[1] = %.16e\n", scoreYmXB[n]);*/
	// residual sigma
	/*double *scoreYmXB_tmp = (double *) fgls_malloc ( n * y_inc * sizeof(double) );*/
	/*dgemm_(TRANS, NO_TRANS,*/
	/*&n, &y_inc, &n, */
	/*&ONE, loops_t->Z, &n, scoreYmXB, &n,*/
	/*&ZERO, scoreYmXB_tmp, &n);*/
	/*for ( ll = 0; ll < y_inc; ll++ )*/
	/*{*/
	/*for ( k = 0; k < n; k++ )*/
	/*{*/
	/*scoreYmXB_tmp[ ll * n + k ] = scoreYmXB_tmp[ ll * n + k ] * Winv[ ll * n + k ];*/
	/*}*/
	/*res_sigma[ll] = ddot_(&n, &scoreYmXB[ ll * n ], &iONE, */
	/*&scoreYmXB[ ll * n ], &iONE) / (n - wXL);*/
	/*printf("sigma[%d] = %.16e\n", ll, res_sigma[ll]);*/
	/*}*/

#if VAMPIR
    VT_USER_END("COMP_LOOP_Y_CODE");
#endif
    for (ib = 0; ib < m; ib += x_b) 
	{
		/*printf("Iter ib: %d\n", ib);*/
		/*printf("%d Iter: %d\n", id, iter);*/
#if VAMPIR
		VT_USER_START("READ_X");
#endif
		/* Read the next X */
	    struct aiocb *aiocb_x_next_p = (struct aiocb *)aiocb_x_next_l[0];
		fgls_aio_read( aiocb_x_next_p,
				       fileno( loops_t->XR_fp ), x_next,
					   (ib + x_b) >= m ? (size_t)MIN( x_b, m ) * wXR * n * sizeof(double) : MIN( (size_t)x_b * wXR * n, (size_t)(m - (ib + x_b)) * wXR * n ) * sizeof(double),
					   (ib + x_b) >= m ? 0 : (off_t)(ib + x_b) * wXR * n * sizeof(double) );

#if VAMPIR
		VT_USER_END("READ_X");
#endif

#if VAMPIR
		VT_USER_START("WAIT_X");
#endif
		/* Wait until the current X is available */
		fgls_aio_suspend( aiocb_x_cur_l, 1, NULL );
#if VAMPIR
		VT_USER_END("WAIT_X");
#endif
#if VAMPIR
      VT_USER_START("COMP_LOOP_X_CODE");
#endif
      x_inc = MIN(x_b, m - ib);
	  for ( j = 0; j < y_inc; j++ )
	  {
		  /*printf("Iter j: %d\n", j);*/
		  /* XR := sqrt(Winv) * Z' * XR */
		  for ( k = 0; k < x_inc*wXR; k++ )
			  for ( l = 0; l < n; l++ )
				  x_copy[ k * n + l ] = x_cur[ k * n + l ] * Winv[j*n + l];
		  
		  /* B_b := WinvZtX * y */
		  /*dgemv_("T", &n, &mpb_real, &ONE, x_cur, &n, y_cur, &iONE, &ZERO, b_cur, &iONE);*/
		  
		  for (i = 0; i < x_inc; i++)
		  {
			  Bij = &b_cur[ j*p*x_inc   + i*p];
			  Vij = &v_cur[ j*p*p*x_inc + i*p*p];

			  /* Building B */
			  // Copy B_T
			  memcpy(Bij, &B_t[j * wXL], wXL * sizeof(double));
			  // B_B := XR' * y
			  dgemv_("T", 
					  &n, &wXR, 
					  &ONE, &x_copy[i * wXR * n], &n, &y_cur[j* n], &iONE, 
					  &ZERO, &Bij[wXL], &iONE);

				/* Building V */
			    // Copy V_TL
				for( k = 0; k < wXL; k++ )
					dcopy_(&wXL, &V_tl[j*wXL*wXL + k*wXL], &iONE, &Vij[k*p], &iONE); // TL
				// V_BL := XR' * XL
				dgemm_("T", "N",
						&wXR, &wXL, &n, // m, n, k
						&ONE, &x_copy[i * wXR * n], &n, &XL[j*wXL*n], &n, // KR KL'
						&ZERO, &Vij[wXL], &p); // BL
				// V_BR := XR' * XR
				dsyrk_("L", "T", 
						&wXR, &n, 
						&ONE, &x_copy[i * wXR * n], &n, 
						&ZERO, &Vij[wXL * p + wXL], &p);

				/* B := inv(V) * y */
				dpotrf_(LOWER, &p, Vij, &p, &info);
				if (info != 0)
				{
					char err[STR_BUFFER_SIZE];
					snprintf(err, STR_BUFFER_SIZE, "dpotrf failed (info: %d)", info);
					error_msg(err, 1);
				}
				dtrsv_(LOWER, NO_TRANS, NON_UNIT, &p, Vij, &p, Bij, &iONE);
				dtrsv_(LOWER,    TRANS, NON_UNIT, &p, Vij, &p, Bij, &iONE);

				// V = res_sigma * inv( X' inv(M) X)
				dpotri_(LOWER, &p, Vij, &p, &info);
				if (info != 0)
				{
					char err[STR_BUFFER_SIZE];
					snprintf(err, STR_BUFFER_SIZE, "dpotri failed (info: %d)", info);
					error_msg(err, 1);
				}
				int p2 = p*p;
				dscal_(&p2, &res_sigma[jb+j], Vij, &iONE);
				for ( k = 0; k < p; k++ )
					Vij[k*p+k] = sqrt(Vij[k*p+k]);
		  }
	  }
#if VAMPIR
      VT_USER_END("COMP_LOOP_X_CODE");
#endif
		
#if VAMPIR
      VT_USER_START("WAIT_B");
#endif
	  /* Wait until previously computed B and V are written */
	  if ( iter > 0)
	  {
		  /*if ( aio_suspend( aiocb_b_prev_l, prev_y_inc, NULL ) != 0 ) // FIX? */
		  /*if ( aio_suspend( aiocb_b_prev_l, y_b, NULL ) != 0 ) // FIX? */
		  for ( k = 0; k < y_b; k++ )
		  {
			fgls_aio_suspend( &aiocb_b_prev_l[k], 1, NULL );
			fgls_aio_suspend( &aiocb_v_prev_l[k], 1, NULL );
		  }
	  }
#if VAMPIR
      VT_USER_END("WAIT_B");
#endif

#if VAMPIR
      VT_USER_START("WRITE_B");
#endif
	  /* Write current B and V */
	  struct aiocb *aiocb_b_cur_p;
	  struct aiocb *aiocb_v_cur_p;
	  for ( k = 0; k < y_inc; k++ )
	  {
		  aiocb_b_cur_p = (struct aiocb *) aiocb_b_cur_l[k];
		  fgls_aio_write( aiocb_b_cur_p,
				          fileno( loops_t->B_fp ), &b_cur[ k * x_inc * p],
						  (size_t)x_inc * p * sizeof(double),
						  ((off_t)(jb+k) * m * p + ib * p) * sizeof(double) );

		  aiocb_v_cur_p = (struct aiocb *) aiocb_v_cur_l[k];
		  fgls_aio_write( aiocb_v_cur_p,
				          fileno( loops_t->V_fp ), &v_cur[ k * x_inc * p * p],
						  (size_t)x_inc * p * p * sizeof(double),
						  ((off_t)(jb+k) * m * p * p + ib * p * p) * sizeof(double) );
	  }
#if VAMPIR
      VT_USER_END("WRITE_B");
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
    gettimeofday(&t1, NULL);
	/*if (id == 0)*/
	/*{*/
	/*printf("Iter (%d - %d) time: %ld ms\n", jb, jb+y_inc, get_diff_ms(&t0, &t1));*/
	/*fflush(stdout);*/
	/*}*/
	/* Swap buffers */
	swap_aiocb( &aiocb_y_cur_l, &aiocb_y_next_l );
	swap_buffers( &y_cur, &y_next);
  }
	/* Wait for the remaining IO operations issued */
#if VAMPIR
      VT_USER_START("WAIT_X");
#endif
  fgls_aio_suspend( aiocb_x_cur_l, 1, NULL );
#if VAMPIR
      VT_USER_END("WAIT_X");
#endif
#if VAMPIR
      VT_USER_START("WAIT_Y");
#endif
  fgls_aio_suspend( aiocb_y_cur_l, 1, NULL );
#if VAMPIR
      VT_USER_END("WAIT_Y");
#endif
#if VAMPIR
      VT_USER_START("WAIT_B");
#endif
  for ( i = 0; i < y_inc; i++ )
  {
    fgls_aio_suspend( &aiocb_b_prev_l[i], 1, NULL );
    fgls_aio_suspend( &aiocb_v_prev_l[i], 1, NULL );
  }
#if VAMPIR
      VT_USER_END("WAIT_B");
#endif

	/* Clean-up */
	  free( x_copy );

  /*free( scoreB_t  );*/
  /*free( scoreV_tl );*/
  /*free( scoreYmXB );*/
  /*free( res_sigma );*/
  
  free( aiocb_b_prev );
  free( aiocb_b_cur  );
  free( aiocb_v_prev );
  free( aiocb_v_cur  );
  free( aiocb_x_cur_l  );
  free( aiocb_x_next_l );
  free( aiocb_y_cur_l  );
  free( aiocb_y_next_l );
  free( aiocb_b_prev_l );
  free( aiocb_b_cur_l  );
  free( aiocb_v_prev_l );
  free( aiocb_v_cur_l  );

  pthread_exit(NULL);
}

/*
 * Performs the computation outside the loops:
 *   - Z W Z' = Phi
 *   - Z' X
 *   - Z' Y
 *   - residual sigma's
 */
int preloop(FGLS_config_t *cf, double *Phi, double *Z, double *W, double *res_sigma) 
{
	/* Threads for ooc gemms */
	ooc_gemm_t gemm_t;
	pthread_t compute_thread;
	int iret;

	/* Buffer sizes for ooc gemms */
	long int chunk_size = 1L << 26; // 64MElems - 28; // 256 MElems
	chunk_size = chunk_size - chunk_size % (cf->n * sizeof(double));
	int num_cols = chunk_size / (cf->n * sizeof(double));

	/* residual sigma's */
	FILE *fp;
	ooc_res_sigma_t res_sigma_t;

	/* Z W Z' = Phi */
	/*printf("\nEigendecomposition of Phi...");*/
	/*fflush(stdout);*/

	eigenDec( cf->n, Phi, Z, W );
	
	/*printf(" Done\n");*/
	/*fflush(stdout);*/

	/* OOC gemms */
	gemm_t.in[0]  = ( double * ) fgls_malloc ( chunk_size * sizeof(double) );
	gemm_t.in[1]  = ( double * ) fgls_malloc ( chunk_size * sizeof(double) );
	gemm_t.out[0] = ( double * ) fgls_malloc ( chunk_size * sizeof(double) );
	gemm_t.out[1] = ( double * ) fgls_malloc ( chunk_size * sizeof(double) );
	/*if ( gemm_t.out[1] == NULL )*/
	/*{*/
	/*fprintf(stderr, "Not enough memory for OOC gemm's\n");*/
	/*exit(EXIT_FAILURE);*/
	/*}*/

	/*printf("Computing Z' XL...");*/
	/*fflush(stdout);*/

	gemm_t.Z = Z;

	/* Set up dimensions */
	gemm_t.m = cf->n;
	gemm_t.n = cf->wXL;
	gemm_t.k = cf->n;
	gemm_t.n_cols_per_buff = num_cols;
	gemm_t.fp_in  = fopen( cf->XL_path, "rb" );
	gemm_t.fp_out = fopen( cf->ZtXL_path, "wb" );
	/*fd_in	= open( cf->X_path, O_RDONLY | O_SYNC );*/
	/*fd_out = open( cf->ZtX_path, O_WRONLY | O_CREAT | O_SYNC );*/
	/*gemm_t.fp_in	= fdopen( fd_in,	"rb" );*/
	/*gemm_t.fp_out = fdopen( fd_out, "wb" );*/

	iret = pthread_create(&compute_thread, NULL, ooc_gemm, (void*)&gemm_t);
	if (iret)
	{
		fprintf(stderr, __FILE__ ": Error creating Computation thread (1): %d\n", iret);
		exit(EXIT_FAILURE);
	}

	pthread_join(compute_thread, NULL);

	fclose( gemm_t.fp_in );
	fclose( gemm_t.fp_out );

	/*printf(" Done\n");*/
	/*fflush(stdout);*/

	/*printf("Computing Z' XR...");*/
	/*fflush(stdout);*/

	/* Set up dimensions */
	gemm_t.m = cf->n;
	gemm_t.n = cf->m * cf->wXR;
	gemm_t.k = cf->n;
	gemm_t.n_cols_per_buff = num_cols;
	gemm_t.fp_in  = fopen( cf->XR_path, "rb" );
	gemm_t.fp_out = fopen( cf->ZtXR_path, "wb" );

	iret = pthread_create(&compute_thread, NULL, ooc_gemm, (void*)&gemm_t);
	if (iret)
	{
		fprintf(stderr, __FILE__ ": Error creating Computation thread (2): %d\n", iret);
		exit(EXIT_FAILURE);
	}

	pthread_join(compute_thread, NULL);

	fclose( gemm_t.fp_in );
	fclose( gemm_t.fp_out );

	/*printf(" Done\n");*/
	/*fflush(stdout);*/

	/*printf("Computing Z' Y...");*/
	/*fflush(stdout);*/

	/* Set up dimensions */
	gemm_t.m = cf->n;
	gemm_t.n = cf->t;
	gemm_t.k = cf->n;
	gemm_t.n_cols_per_buff = num_cols;
	gemm_t.fp_in  = fopen( cf->Y_path, "rb" );
	gemm_t.fp_out = fopen( cf->ZtY_path, "wb" );

	iret = pthread_create(&compute_thread, NULL, ooc_gemm, (void*)&gemm_t);
	if (iret)
	{
		fprintf(stderr, __FILE__ ": Error creating Computation thread (3): %d\n", iret);
		exit(EXIT_FAILURE);
	}

	pthread_join(compute_thread, NULL);

	fclose( gemm_t.fp_in );
	fclose( gemm_t.fp_out );

	/*printf(" Done\n");*/
	/*fflush(stdout);*/

	/* Clean-up */
	free(gemm_t.in[0]);
	free(gemm_t.in[1]);
	free(gemm_t.out[0]);
	free(gemm_t.out[1]);

	/* 
	 * residual sigma's 
	 */
	res_sigma_t.h2      = (double*) fgls_malloc ( cf->t * sizeof(double) );
	res_sigma_t.sigma2  = (double*) fgls_malloc ( cf->t * sizeof(double) );
	res_sigma_t.XL_orig = (double*) fgls_malloc ( cf->wXL * cf->n * sizeof(double) );
	res_sigma_t.ZtXL    = (double*) fgls_malloc ( cf->wXL * cf->n * sizeof(double) );

	res_sigma_t.Z = Z;
	res_sigma_t.W = W;
	res_sigma_t.fp_Y   = fopen( cf->Y_path, "rb" );
	res_sigma_t.fp_ZtY = fopen( cf->ZtY_path, "rb" );
	res_sigma_t.res_sigma = res_sigma;
	/* Load h2 */
	fp = fopen( cf->h_path, "rb" );
	sync_read( res_sigma_t.h2, fp, cf->t, 0 );
	fclose( fp );
	/* Load sigma2 */
	fp = fopen( cf->sigma_path, "rb" );
	sync_read( res_sigma_t.sigma2, fp, cf->t, 0 );
	fclose( fp );
	/* Load XL original */
	fp = fopen( cf->XL_path, "rb" );
	sync_read( res_sigma_t.XL_orig, fp, cf->n * cf->wXL, 0 );
	fclose( fp );
	/* Load ZtXl */
	fp = fopen( cf->ZtXL_path, "rb" );
	sync_read( res_sigma_t.ZtXL, fp, cf->n * cf->wXL, 0 );
	fclose( fp );

	residual_sigma( cf, &res_sigma_t );

	/* Clean-up */
	free( res_sigma_t.h2 );
	free( res_sigma_t.sigma2 );
	free( res_sigma_t.XL_orig );
	free( res_sigma_t.ZtXL );

	fclose( res_sigma_t.fp_Y );
	fclose( res_sigma_t.fp_ZtY );

	return 0;
}

/* 
 * Eigendecomposition of Phi
 *
 * Z W Z' = Phi
 */
void eigenDec(int n, double *Phi, double *Z, double *W)
{
	int nb = 192;
	int idummy, nCompPairs, *isuppz, *iwork, info = 0,
	    lwork = n * (nb + 6),
	    liwork = 10 * n;
	double ddummy = -1.0, *work;

	work   = (double *) fgls_malloc ( lwork * sizeof(double) );
	iwork  = (int *)    fgls_malloc ( liwork * sizeof(int) );
	isuppz = (int *)    fgls_malloc ( 2 * n * sizeof(int) );

	dsyevr_("V", "A", "L", &n, Phi, &n, 
	        &ddummy, &ddummy, &idummy, &idummy, &ddummy, 
	        &nCompPairs, W, Z, &n, isuppz, 
	        work, &lwork, iwork, &liwork, &info);

	free( work );
	free( iwork );
	free( isuppz );

	if (info != 0)
	{
		char err[STR_BUFFER_SIZE];
		snprintf(err, STR_BUFFER_SIZE, "dsyevr failed (info: %d)", info);
		error_msg(err, 1);
	}
}

/*
 * Out-of-core gemms:
 *   - Z' XL
 *   - Z' XR
 *   - Z' Y
 */
void* ooc_gemm( void *in ) 
{
	ooc_gemm_t *gemm_t = ( ooc_gemm_t *)in;

	/* Problem dimensions */
	int m = gemm_t->m,
	    n = gemm_t->n,
	    k = gemm_t->k,
	    n_cols_per_buff = gemm_t->n_cols_per_buff;
	long int max_elems = n_cols_per_buff * m;

	/* Double buffering */
	double *in_cur   = gemm_t->in[0];
	double *in_next  = gemm_t->in[1];
	double *out_prev = gemm_t->out[0];
	double *out_cur  = gemm_t->out[1];

	/* BLAS constants */
	double ONE  = 1.0;
	double ZERO = 0.0;

	/* Asynchronous IO data structures */
	struct aiocb aiocb_in_cur,   aiocb_in_next,
	             aiocb_out_prev, aiocb_out_cur;

	const struct aiocb ** aiocb_in_cur_l,
	                   ** aiocb_in_next_l,
	                   ** aiocb_out_prev_l,
	                   ** aiocb_out_cur_l;

	aiocb_in_cur_l   = (const struct aiocb **) fgls_malloc (sizeof(struct aiocb *));
	aiocb_in_next_l  = (const struct aiocb **) fgls_malloc (sizeof(struct aiocb *));
	aiocb_out_prev_l = (const struct aiocb **) fgls_malloc (sizeof(struct aiocb *));
	aiocb_out_cur_l  = (const struct aiocb **) fgls_malloc (sizeof(struct aiocb *));

	aiocb_in_cur_l[0]   = &aiocb_in_cur;
	aiocb_in_next_l[0]  = &aiocb_in_next;
	aiocb_out_prev_l[0] = &aiocb_out_prev;
	aiocb_out_cur_l[0]  = &aiocb_out_cur;

	/* Read first piece of "in" */
	fgls_aio_read( &aiocb_in_cur,
	               fileno( gemm_t->fp_in ), in_cur,
	               MIN( (size_t)max_elems, (size_t)k * n ) * sizeof(double), 0);

	int cur_n;
	int i;
	for ( i = 0; i < n; i += n_cols_per_buff ) 
	{
		/* Read next piece of "in" */
		struct aiocb *aiocb_in_next_p = (struct aiocb *)aiocb_in_next_l[0];
		fgls_aio_read( aiocb_in_next_p,
					   fileno( gemm_t->fp_in ), in_next,
					   i + n_cols_per_buff > n ? 0 : MIN( max_elems, ( n - (size_t)( i + n_cols_per_buff ) ) * k ) * sizeof(double),
					   i + n_cols_per_buff > n ? 0 : (off_t)(i + n_cols_per_buff) * k * sizeof(double) );

		/* Wait for current piece of "in" */
		fgls_aio_suspend( aiocb_in_cur_l, 1, NULL );

		/* Compute */
		cur_n = MIN( n_cols_per_buff, (n - i) );
		dgemm_("T", "N", &m, &cur_n, &k, &ONE, gemm_t->Z, &m, in_cur, &m, &ZERO, out_cur, &m);

		/* Wait until previous piece of "out" is written */
		if ( i > 0)
		{
			fgls_aio_suspend( aiocb_out_prev_l, 1, NULL );
		}
		/* Write current piece of "out" */
		struct aiocb *aiocb_out_cur_p = (struct aiocb *)aiocb_out_cur_l[0];
		fgls_aio_write( aiocb_out_cur_p,
						fileno( gemm_t->fp_out ), out_cur,
						MIN( max_elems, (size_t)(n - i) * m ) * sizeof(double),
						(off_t)i * m * sizeof(double) );

		/* Swap buffers */
		swap_aiocb( &aiocb_in_cur_l,   &aiocb_in_next_l );
		swap_aiocb( &aiocb_out_prev_l, &aiocb_out_cur_l );
		swap_buffers( &in_cur,   &in_next );
		swap_buffers( &out_prev, &out_cur );
	}

	/* Wait for the remaining io calls issued */
	fgls_aio_suspend( aiocb_in_cur_l, 1, NULL );
	fgls_aio_suspend( aiocb_out_prev_l, 1, NULL );

	free( aiocb_in_cur_l   );
	free( aiocb_in_next_l  );
	free( aiocb_out_prev_l );
	free( aiocb_out_cur_l  );

	pthread_exit(NULL);
}

/*
 * Compute residual sigma's
 */
void residual_sigma( FGLS_config_t *fgls_cf, ooc_res_sigma_t *cf ) 
{
	/* Problem dimensions */
	int n = fgls_cf->n,
	    t = fgls_cf->t,
		wXL = fgls_cf->wXL,
	    y_b = fgls_cf->y_b,
		j, k, l, ll;

	/* BLAS constants */
	double ZERO = 0.0;
	double ONE  = 1.0;
	double MINUS_ONE  = -1.0;
	int iONE = 1;
	/* LAPACK error */
	int info;

	/* intermediate values */
	double *alpha = (double *) fgls_malloc (y_b * sizeof(double));
	double *beta  = (double *) fgls_malloc (y_b * sizeof(double));
	double *Winv  = (double *) fgls_malloc (n * y_b * sizeof(double));
	double *XL_copy = (double *) fgls_malloc (n * wXL * y_b * sizeof(double));
	/*memcpy( XL_copy, cf->XL_orig, wXL * n * sizeof(double) );*/

	/* sigma2.score */
	double *scoreB_t  = (double *) fgls_malloc (wXL * y_b * sizeof(double));
	double *scoreV_tl = (double *) fgls_malloc (wXL * wXL * y_b * sizeof(double));
	double *scoreYmXB = (double *) fgls_malloc (n * y_b * sizeof(double));
	double *scoreYmXB_tmp = (double *) fgls_malloc ( n * y_b * sizeof(double) );
	/*double *res_sigma = (double *) fgls_malloc (y_b * sizeof(double));*/

	/* Double buffering */
	double *Y_cur    = (double *) fgls_malloc (n * y_b * sizeof(double));
	double *Y_next   = (double *) fgls_malloc (n * y_b * sizeof(double));
	double *ZtY_cur  = (double *) fgls_malloc (n * y_b * sizeof(double));
	double *ZtY_next = (double *) fgls_malloc (n * y_b * sizeof(double));

	/* Asynchronous IO data structures */
	struct aiocb aiocb_Y_cur,   aiocb_Y_next,
	             aiocb_ZtY_cur, aiocb_ZtY_next;

	const struct aiocb ** aiocb_Y_cur_l,
	                   ** aiocb_Y_next_l,
	                   ** aiocb_ZtY_cur_l,
	                   ** aiocb_ZtY_next_l;

	aiocb_Y_cur_l   = (const struct aiocb **) fgls_malloc (sizeof(struct aiocb *));
	aiocb_Y_next_l  = (const struct aiocb **) fgls_malloc (sizeof(struct aiocb *));
	aiocb_ZtY_cur_l = (const struct aiocb **) fgls_malloc (sizeof(struct aiocb *));
	aiocb_ZtY_next_l  = (const struct aiocb **) fgls_malloc (sizeof(struct aiocb *));

	aiocb_Y_cur_l[0]    = &aiocb_Y_cur;
	aiocb_Y_next_l[0]   = &aiocb_Y_next;
	aiocb_ZtY_cur_l[0]  = &aiocb_ZtY_cur;
	aiocb_ZtY_next_l[0] = &aiocb_ZtY_next;

	/* Read first piece of "Y" */
	fgls_aio_read( &aiocb_Y_cur,
	               fileno( cf->fp_Y ), Y_cur,
	               MIN( (size_t)y_b, (size_t)t) * n * sizeof(double), 0);
	/* Read first piece of "ZtY" */
	fgls_aio_read( &aiocb_ZtY_cur,
	               fileno( cf->fp_ZtY ), ZtY_cur,
	               MIN( (size_t)y_b, (size_t)t) * n * sizeof(double), 0);
	int y_inc;
	for ( j = 0; j < t; j+=y_b )
	{
		y_inc = MIN( y_b, t - j );
		/* Read next piece of Y */
		struct aiocb *aiocb_Y_next_p = (struct aiocb *)aiocb_Y_next_l[0];
		fgls_aio_read( aiocb_Y_next_p,
					   fileno( cf->fp_Y ), Y_next,
					   j + y_b > t ? 0 : MIN( y_b, t - (size_t)( j + y_b ) ) * n * sizeof(double),
					   j + y_b > t ? 0 : (off_t)(j + y_b) * n * sizeof(double) );
		/* Read next piece of ZtY */
		struct aiocb *aiocb_ZtY_next_p = (struct aiocb *)aiocb_ZtY_next_l[0];
		fgls_aio_read( aiocb_ZtY_next_p,
					   fileno( cf->fp_ZtY ), ZtY_next,
					   j + y_b > t ? 0 : MIN( y_b, t - (size_t)( j + y_b ) ) * n * sizeof(double),
					   j + y_b > t ? 0 : (off_t)(j + y_b) * n * sizeof(double) );

		/* Wait for current piece of Y */
		fgls_aio_suspend( aiocb_Y_cur_l, 1, NULL );
		/* Wait for current piece of Y */
		fgls_aio_suspend( aiocb_ZtY_cur_l, 1, NULL );

		/* Compute */
		for (k = 0; k < y_inc; k++)
		{
			alpha[k] = cf->sigma2[j + k] * cf->h2[j + k];
			beta[k]  = cf->sigma2[j + k] * (1 - cf->h2[j + k]);
		}

		/* Winv := sqrt( inv( alpha W - beta I ) ) */
		// Best order? sqrt - inv
		//
		// Possibly GER
		for (k = 0; k < y_inc; k++)
			for (l = 0; l < n; l++)
				Winv[k*n + l] = sqrt(1.0 / (alpha[k] * cf->W[l] + beta[k]));

		/* y := sqrt(Winv) * Z' * y */
		for (k = 0; k < y_inc; k++)
			for (l = 0; l < n; l++)
				ZtY_cur[k*n+l] *= Winv[k*n+l];

		/* XL := sqrt(Winv) * Z' * XL */
		for (ll = 0; ll < y_inc; ll++)
		  for ( k = 0; k < wXL; k++ )
			  for ( l = 0; l < n; l++ )
				  XL_copy[ ll * wXL * n + k * n + l ] = cf->ZtXL[ k * n + l ] * Winv[ll * n + l];
		  
		/* B_t := XL' * y */
		for (ll = 0; ll < y_inc; ll++)
		  dgemv_("T", &n, &wXL, &ONE, &XL_copy[ll * wXL * n], &n, &ZtY_cur[ll * n], &iONE, &ZERO, &scoreB_t[ll * wXL], &iONE);

		/* V_tl := Compute Top-left part of V */
		for (ll = 0; ll < y_inc; ll++)
		{
			dsyrk_("L", "T", // LOWER, NO_TRANS, 
					&wXL, &n, // n, k
					&ONE, &XL_copy[ll * wXL * n], &n, // KL KL' 
					&ZERO, &scoreV_tl[ll * wXL * wXL], &wXL); // V_TL
		}

		/* Compute sigma2.score's */
		// bt
		double *sc_B_t, *sc_V_tl;
		for ( ll = 0; ll < y_inc; ll++ )
		{
			sc_B_t  = &scoreB_t [ ll * wXL ];
			sc_V_tl = &scoreV_tl[ ll * wXL * wXL ];

			dpotrf_(LOWER, &wXL, sc_V_tl, &wXL, &info);
			if (info != 0)
			{
				char err[STR_BUFFER_SIZE];
				snprintf(err, STR_BUFFER_SIZE, "sigma2.score: dpotrf(scoreV) failed (info: %d) - j: %d", info, j+ll);
				error_msg(err, 1);
			}
			dtrsv_(LOWER, NO_TRANS, NON_UNIT, &wXL, sc_V_tl, &wXL, sc_B_t, &iONE);
			dtrsv_(LOWER,    TRANS, NON_UNIT, &wXL, sc_V_tl, &wXL, sc_B_t, &iONE);
		}
		// YmXB
		memcpy( scoreYmXB, Y_cur, n * y_inc * sizeof(double) );
		dgemm_(NO_TRANS, NO_TRANS,
				&n, &y_inc, &wXL, 
				&MINUS_ONE, cf->XL_orig, &n, scoreB_t, &wXL,
				&ONE, scoreYmXB, &n);
		// residual sigma
		dgemm_(TRANS, NO_TRANS,
				&n, &y_inc, &n, 
				&ONE, cf->Z, &n, scoreYmXB, &n,
				&ZERO, scoreYmXB_tmp, &n);
		for ( ll = 0; ll < y_inc; ll++ )
		{
			for ( k = 0; k < n; k++ )
			{
				scoreYmXB_tmp[ ll * n + k ] = scoreYmXB_tmp[ ll * n + k ] * Winv[ ll * n + k ];
			}
			cf->res_sigma[j + ll] = ddot_(&n, &scoreYmXB_tmp[ ll * n ], &iONE, 
									  		  &scoreYmXB_tmp[ ll * n ], &iONE) / (n - wXL);
		}

		/* Swap buffers */
		swap_aiocb( &aiocb_Y_cur_l,   &aiocb_Y_next_l );
		swap_aiocb( &aiocb_ZtY_cur_l, &aiocb_ZtY_next_l );
		swap_buffers( &Y_cur,   &Y_next );
		swap_buffers( &ZtY_cur, &ZtY_next );
	}

	/* Wait for the remaining io calls issued */
	fgls_aio_suspend( aiocb_Y_cur_l, 1, NULL );
	fgls_aio_suspend( aiocb_ZtY_cur_l, 1, NULL );

	/* Clean-up */
	free( alpha );
	free( beta );
	free( Winv );
	free( XL_copy );

	free( scoreB_t );
	free( scoreV_tl );
	free( scoreYmXB );
	free( scoreYmXB_tmp );

	free( Y_cur );
	free( Y_next );
	free( ZtY_cur );
	free( ZtY_next );

	free( aiocb_Y_cur_l    );
	free( aiocb_Y_next_l   );
	free( aiocb_ZtY_cur_l  );
	free( aiocb_ZtY_next_l );
	
	return;
}
