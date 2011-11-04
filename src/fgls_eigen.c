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
#include "fgls_eigen.h"

#if VAMPIR
  #include "vt_user.h"
#endif

void eigenDec(int n, double *Phi, double *Z, double *W);
void* compute_thread_func(void* in);
void* ooc_gemm_comp( void *in );
/*void* ooc_gemm_io( void *in );*/

/*
 * Entry point for the solution of the Feasible Generalized Least-Squares problem
 */
int fgls_eigen(int n, int p, int m, int t, int wXL, int wXR,
               int x_b, int y_b, int num_threads,
			   char *Phi_path, char *h2_path, char *sigma2_path,
			   char *XL_path, char *XR_path, char *Y_path,
			   char *B_path, char *V_path)
{
	FGLS_config_t cf;

  int rc, i;
  char numths_str[STR_BUFFER_SIZE];
  pthread_t *compute_threads;

  double *Phi,
		 *Z;
  FILE *Phi_fp, *h_fp, *sigma_fp,
	   *XL_fp;
  ooc_loops_t loops_t;
  ooc_loops_t *loops_t_comp;

  /* Check input values */
  /*printf("n: %d\np: %d\nm: %d\nt: %d\nwXL: %d\nwXR: %d\n", n, p, m, t, wXL, wXR);*/
  /*printf("x_b: %d\ny_b: %d\nnths: %d\n", x_b, y_b, num_threads);*/
  /*printf("Phi: %s\nh2: %s\ns2: %s\n", cf.Phi_path, cf.h_path, cf.sigma_path);*/
  /*printf("XL: %s\nXR: %s\nY: %s\n", cf.XL_path, cf.XR_path, cf.Y_path);*/
  /*printf("B: %s\nV: %s\n", cf.B_path, cf.V_path);*/

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
  Z   = ( double * ) fgls_malloc ( cf.n * cf.n * sizeof(double) );
  loops_t.W = ( double * ) fgls_malloc ( cf.n * sizeof(double) );
  /* Load Phi */
  Phi_fp = fopen( cf.Phi_path, "rb" );
  sync_read( Phi, Phi_fp, cf.n * cf.n, 0 );
  fclose( Phi_fp );
  /* Compute the pre-loop operations */
  snprintf(numths_str, STR_BUFFER_SIZE, "%d", cf.NUM_COMPUTE_THREADS);
  setenv("OMP_NUM_THREADS", numths_str, 1);
  preloop(&cf, Phi, Z, loops_t.W);
  setenv("OMP_NUM_THREADS", "1", 1);
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
  // Double-buffering for the operands that do not fit in RAM
  for (i = 0; i < NUM_BUFFERS_PER_THREAD; i++) 
  {
    loops_t.X[i] = ( double* ) fgls_malloc ( cf.NUM_COMPUTE_THREADS * cf.x_b * cf.p * cf.n * sizeof(double) );
    loops_t.Y[i] = ( double* ) fgls_malloc ( cf.NUM_COMPUTE_THREADS * cf.y_b * cf.n * sizeof(double) );
    loops_t.B[i] = ( double* ) fgls_malloc ( cf.NUM_COMPUTE_THREADS * cf.x_b * cf.y_b * cf.p * sizeof(double) );
    loops_t.V[i] = ( double* ) fgls_malloc ( cf.NUM_COMPUTE_THREADS * cf.x_b * cf.y_b * cf.p * cf.p * sizeof(double) );
  }
  // The original in XL[0]
  // The copy to overwrite in XL[1]
  loops_t.XL[0] = ( double* ) fgls_malloc ( cf.wXL * cf.n * sizeof(double) );
  loops_t.XL[1] = ( double* ) fgls_malloc ( cf.NUM_COMPUTE_THREADS * cf.wXL * cf.n * cf.y_b * sizeof(double) );
  // Reusable b_Top and V_TopLeft
  loops_t.XL_b  = ( double* ) fgls_malloc ( cf.NUM_COMPUTE_THREADS * cf.wXL * cf.y_b * sizeof(double) );
  loops_t.XLtXL = ( double* ) fgls_malloc ( cf.NUM_COMPUTE_THREADS * cf.wXL * cf.wXL * cf.y_b * sizeof(double) );

  loops_t.h     = ( double* ) fgls_malloc ( cf.t * sizeof(double) );
  loops_t.sigma = ( double* ) fgls_malloc ( cf.t * sizeof(double) );

  loops_t.alpha   = ( double* ) fgls_malloc ( cf.NUM_COMPUTE_THREADS * cf.y_b * sizeof(double) );
  loops_t.beta    = ( double* ) fgls_malloc ( cf.NUM_COMPUTE_THREADS * cf.y_b * sizeof(double) );
  loops_t.Winv    = ( double* ) fgls_malloc ( cf.NUM_COMPUTE_THREADS * cf.x_b * cf.n * sizeof(double) );
  /*loops_t.xtSx    = ( double* ) malloc ( cf.NUM_COMPUTE_THREADS * cf.p * cf.p * sizeof(double) );*/

  // Files for XR, Y, B and V
  loops_t.XR_fp = fopen( cf.ZtXR_path, "rb");
  loops_t.Y_fp  = fopen( cf.ZtY_path, "rb");
  loops_t.B_fp  = fopen( cf.B_path, "wb");
  loops_t.V_fp  = fopen( cf.V_path, "wb");

  // Load XL, h and sigma
  XL_fp = fopen( cf.ZtXL_path, "rb" );
  sync_read( loops_t.XL[0], XL_fp, cf.wXL * cf.n, 0 );
  fclose( XL_fp );
  h_fp = fopen( cf.h_path, "r");
  sync_read(loops_t.h, h_fp, cf.t, 0);
  fclose( h_fp );
  sigma_fp = fopen( cf.sigma_path, "r");
  sync_read(loops_t.sigma, sigma_fp, cf.t, 0);
  fclose( sigma_fp );

  printf("LOOPS\n");

  // Spawn threads
  loops_t_comp = ( ooc_loops_t* ) fgls_malloc ( cf.NUM_COMPUTE_THREADS * sizeof(ooc_loops_t) );
  compute_threads = ( pthread_t * ) fgls_malloc ( cf.NUM_COMPUTE_THREADS * sizeof(pthread_t) );
  /*if (compute_threads == NULL)*/
  /*{*/
  /*fprintf(stderr, "Not enough memory\n");*/
  /*exit(EXIT_FAILURE);*/
  /*}*/
  for (i = 0; i < cf.NUM_COMPUTE_THREADS; i++) 
  {
    memcpy((void*)&loops_t_comp[i], (void*)&loops_t, sizeof(ooc_loops_t));    
    loops_t_comp[i].id = i;

    rc = pthread_create(&compute_threads[i], NULL, compute_thread_func, (void*)&loops_t_comp[i]);
    if (rc) 
	{
      fprintf(stderr, "ERROR in " __FILE__ ": pthread_create() returned %d\n", rc);
      pthread_exit(NULL);
    }
  }

  void* retval;
  for (i = 0; i < cf.NUM_COMPUTE_THREADS; i++)
    pthread_join(compute_threads[i], &retval);
#if VAMPIR
  VT_USER_START("END");
#endif

  printf("LOOPS Done\n");

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
  free( loops_t.XL[0] );
  free( loops_t.XL[1] );
  free( loops_t.XL_b  );
  free( loops_t.XLtXL );

  free( loops_t.h );
  free( loops_t.sigma );
  free( loops_t.W );

  free( loops_t.alpha );
  free( loops_t.beta );
  free( loops_t.Winv );
  /*free( loops_t.xtSx );*/

  free( loops_t_comp );
  free( compute_threads );

  return 0;
}

void* compute_thread_func(void* in) 
{
  ooc_loops_t *loops_t = ( ooc_loops_t* ) in;
  FGLS_config_t *cf = loops_t->cf;

  int m = cf->m,
	  n = cf->n,
	  p = cf->p,
	  t = cf->t,
	  wXL = cf->wXL,
	  wXR = cf->wXR,
	  x_b = cf->x_b,
	  y_b = cf->y_b,
	  id = loops_t->id;

  double *xl     = &loops_t->XL[1][id * wXL * n * y_b];
  double *xl_b   = &loops_t->XL_b[id * wXL * y_b];
  double *xltxl  = &loops_t->XLtXL[id * wXL * wXL * y_b];
  double *x_cur  = &loops_t->X[0][id * n * wXR * x_b];
  double *x_next = &loops_t->X[1][id * n * wXR * x_b];
  double *y_cur  = &loops_t->Y[0][id * y_b * n];
  double *y_next = &loops_t->Y[1][id * y_b * n];
  double *b_cur  = &loops_t->B[0][id * p * x_b * y_b];
  double *b_prev = &loops_t->B[1][id * p * x_b * y_b];
  double *v_cur  = &loops_t->V[0][id * p * p * x_b * y_b];
  double *v_prev = &loops_t->V[1][id * p * p * x_b * y_b];
  double *h      =  loops_t->h;
  double *sigma  =  loops_t->sigma;
  double *W      =  loops_t->W;

  double *Winv = &loops_t->Winv[id * y_b * n];
  /*double *xtSx = &loops_t->xtSx[id * p * p];*/

  double ONE = 1.0;
  double ZERO = 0.0;
  int iONE = 1;
  /*int iZERO = 0;*/
  double *alpha = loops_t->alpha, 
		 *beta  = loops_t->beta;

  int ib, jb, i, j, k, l, ll;
  int x_inc, y_inc;
  double *Bij, *Vij;
  /*int mpb;*/
  /*int mpb_real;*/
  int info;

  struct aiocb aiocb_x_cur,  aiocb_x_next,
			   aiocb_y_cur,  aiocb_y_next,
			   *aiocb_b_prev, *aiocb_b_cur,
			   *aiocb_v_prev, *aiocb_v_cur;
  double *x_copy = (double *) malloc (x_b * wXR * n * sizeof(double));

  aiocb_b_prev = (struct aiocb *) fgls_malloc (y_b * sizeof(struct aiocb));
  aiocb_b_cur  = (struct aiocb *) fgls_malloc (y_b * sizeof(struct aiocb));
  aiocb_v_prev = (struct aiocb *) fgls_malloc (y_b * sizeof(struct aiocb));
  aiocb_v_cur  = (struct aiocb *) fgls_malloc (y_b * sizeof(struct aiocb));
  struct aiocb **aiocb_x_cur_l,//  = { &aiocb_x_cur }, 
			   **aiocb_x_next_l,// = { &aiocb_x_next },
			   **aiocb_y_cur_l,//  = { &aiocb_y_cur },
			   **aiocb_y_next_l,// = { &aiocb_y_next },
			   **aiocb_b_prev_l,// = {  aiocb_b_prev },
			   **aiocb_b_cur_l,//  = {  aiocb_b_cur };
			   **aiocb_v_prev_l,// = {  aiocb_v_prev },
			   **aiocb_v_cur_l;//  = {  aiocb_v_cur };

  aiocb_x_cur_l  = (struct aiocb **) fgls_malloc (sizeof(struct aiocb *));
  aiocb_x_next_l = (struct aiocb **) fgls_malloc (sizeof(struct aiocb *));
  aiocb_y_cur_l  = (struct aiocb **) fgls_malloc (sizeof(struct aiocb *));
  aiocb_y_next_l = (struct aiocb **) fgls_malloc (sizeof(struct aiocb *));
  aiocb_b_prev_l = (struct aiocb **) fgls_malloc (y_b * sizeof(struct aiocb *));
  aiocb_b_cur_l  = (struct aiocb **) fgls_malloc (y_b * sizeof(struct aiocb *));
  aiocb_v_prev_l = (struct aiocb **) fgls_malloc (y_b * sizeof(struct aiocb *));
  aiocb_v_cur_l  = (struct aiocb **) fgls_malloc (y_b * sizeof(struct aiocb *));

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

  /*bzero( (char *)&aiocb_x_cur,  sizeof(struct aiocb) );*/
  /*bzero( (char *)&aiocb_x_next, sizeof(struct aiocb) );*/
  /*bzero( (char *)&aiocb_y_cur,  sizeof(struct aiocb) );*/
  /*bzero( (char *)&aiocb_y_next, sizeof(struct aiocb) );*/
  /*bzero( (char *)aiocb_b_prev, y_b * sizeof(struct aiocb) );*/
  /*bzero( (char *)aiocb_b_cur,  y_b * sizeof(struct aiocb) );*/
  /*bzero( (char *)aiocb_v_prev, y_b * sizeof(struct aiocb) );*/
  /*bzero( (char *)aiocb_v_cur,  y_b * sizeof(struct aiocb) );*/

  fgls_aio_read( &aiocb_x_cur, 
		         fileno( loops_t->XR_fp ), x_cur,
				 MIN( x_b, m ) * wXR * n * sizeof(double), 0 );
  /*aiocb_x_cur.aio_fildes = fileno( loops_t->XR_fp );*/
  /*aiocb_x_cur.aio_buf = x_cur;*/
  /*aiocb_x_cur.aio_nbytes = MIN( x_b, m ) * wXR * n * sizeof(double);*/
  /*aiocb_x_cur.aio_offset = 0;*/
  /**//*aiocb_x_cur.aio_flag = AIO_RAW;*/
  /*aio_read( &aiocb_x_cur );*/

  fgls_aio_read( &aiocb_y_cur,
		         fileno( loops_t->Y_fp ), y_cur,
				 id >= t ? 0 : n * y_b * sizeof(double),
				 id >= t ? 0 : id * n * y_b * sizeof(double) );
  /*aiocb_y_cur.aio_fildes = fileno( loops_t->Y_fp );*/
  /*aiocb_y_cur.aio_buf = y_cur;*/
  /*aiocb_y_cur.aio_nbytes = id >= t ? 0 : n * y_b * sizeof(double);*/
  /*aiocb_y_cur.aio_offset = id >= t ? 0 : id * n * y_b * sizeof(double);*/
  /**//*aiocb_y_cur.aio_flag = AIO_RAW;*/
  /*aio_read( &aiocb_y_cur );*/

  /*aiocb_x_next.aio_fildes = fileno( loops_t->XR_fp );*/
  /**//*aiocb_x_next.aio_flag = AIO_RAW;*/
  /*aiocb_x_next.aio_buf = x_next;*/
  /*aiocb_y_next.aio_fildes = fileno( loops_t->Y_fp );*/
  /**//*aiocb_y_next.aio_flag = AIO_RAW;*/
  /*aiocb_y_next.aio_buf = y_next;*/

  /*for ( i = 0; i < y_b; i++ )*/
  /*{*/
  /*aiocb_b_prev[i].aio_fildes = fileno( loops_t->B_fp );*/
  /**//*aiocb_b_prev.aio_flag = AIO_RAW;*/
  /*aiocb_b_prev[i].aio_buf = &b_prev[p * x_b * i];*/
  /*aiocb_b_cur[i].aio_fildes = fileno( loops_t->B_fp );*/
  /**//*aiocb_b_cur.aio_flag = AIO_RAW;*/
  /*aiocb_b_cur[i].aio_buf = &b_cur[p * x_b * i];*/

  /*aiocb_v_prev[i].aio_fildes = fileno( loops_t->V_fp );*/
  /**//*aiocb_v_prev.aio_flag = AIO_RAW;*/
  /*aiocb_v_prev[i].aio_buf = &v_prev[p * p * x_b * i];*/
  /*aiocb_v_cur[i].aio_fildes = fileno( loops_t->V_fp );*/
  /**//*aiocb_v_cur.aio_flag = AIO_RAW;*/
  /*aiocb_v_cur[i].aio_buf = &v_cur[p * p * x_b * i];*/
  /*}*/

  int iter = 0;
  struct timeval t0, t1;
  /*int prev_y_inc = 0;*/
  for (jb = id * y_b; jb < t; jb += cf->NUM_COMPUTE_THREADS * y_b) 
  {
	  /*printf("Iter jb: %d\n", jb);*/
    gettimeofday(&t0, NULL);
	y_inc = MIN( y_b, t - jb );
#if VAMPIR
      VT_USER_START("READ_Y");
#endif
    //read( y_cur, loops_t->Y_fp, 
	//		j + cf->NUM_COMPUTE_THREADS > t ? 0 : MIN( cf->NUM_COMPUTE_THREADS * n, (t - (j + cf->NUM_COMPUTE_THREADS)) * n ), 
	//		j + cf->NUM_COMPUTE_THREADS > t ? 0 : (j + cf->NUM_COMPUTE_THREADS) * n );
	struct aiocb *aiocb_y_next_p = (struct aiocb *)aiocb_y_next_l[0];
	fgls_aio_read( aiocb_y_next_p,
			       fileno( loops_t->Y_fp ), y_next,
				   jb + cf->NUM_COMPUTE_THREADS * y_b >= t ? 0 : n * y_b * sizeof(double),
				   jb + cf->NUM_COMPUTE_THREADS * y_b >= t ? 0 : (jb + y_b * cf->NUM_COMPUTE_THREADS) * n * sizeof(double) );

	/*bzero( (char *)aiocb_y_next_p,  sizeof(struct aiocb) );*/
	/**//*printf("Reading Y = %d in %p\n", 1+j+cf->NUM_COMPUTE_THREADS, y_next);*/
	/*aiocb_y_next_p->aio_buf = y_next;*/
	/*aiocb_y_next_p->aio_fildes = fileno( loops_t->Y_fp );*/
	/*aiocb_y_next_p->aio_nbytes = jb + cf->NUM_COMPUTE_THREADS * y_b >= t ? 0 : n * y_b * sizeof(double);*/
	/*aiocb_y_next_p->aio_offset = jb + cf->NUM_COMPUTE_THREADS * y_b >= t ? 0 : (jb + y_b * cf->NUM_COMPUTE_THREADS) * n * sizeof(double);*/
	/**//*aiocb_y.aio_flag = AIO_RAW;*/
	/**//*printf("Nbytes: %d\n",  aiocb_y_next_p->aio_nbytes);*/
	/*aio_read( aiocb_y_next_p );*/
#if VAMPIR
      VT_USER_END("READ_Y");
#endif

#if VAMPIR
      VT_USER_START("WAIT_X");
#endif
	/*printf("xl value: %f\n", *loops_t->XL[0]);*/
	for (ll = 0; ll < y_inc; ll++)
		memcpy( &xl[ll * wXL * n], loops_t->XL[0], wXL * n * sizeof(double) );
#if VAMPIR
      VT_USER_END("WAIT_X");
#endif
#if VAMPIR
    VT_USER_START("WAIT_Y");
#endif
	if ( aio_suspend( aiocb_y_cur_l, 1, NULL ) != 0 )
		perror("Suspend Y error\n");
	/*printf("Return value: %d\n", aio_return((struct aiocb *)aiocb_y_cur_l[0]));*/
#if VAMPIR
    VT_USER_END("WAIT_Y");
#endif
	/*y_cur = (double *)aiocb_y_cur_l[0]->aio_buf;*/
	/*printf("ZtY(1,   %2d)    = %8f\n", j+1, *y_cur);*/
	/*printf("ZtY(500, %2d)    = %8f\n", j+1, y_cur[499]);*/

#if VAMPIR
    VT_USER_START("COMP_LOOP_Y_CODE");
#endif
	printf("h: %f\n", h[0]);
    for (k = 0; k < y_inc; k++)
	{
    	alpha[k] = sigma[jb + k] * h[jb + k]; // * h[j];
    	beta[k]  = sigma[jb + k] * (1 - h[jb + k]); // * h[j]);
	}

    BEGIN_TIMING();
    /* 2) W = sqrt(alpha W - beta I)^-1 */
    // Best order? sqrt - inv
	//
	// Possibly GER
    for (k = 0; k < y_inc; k++)
        for (l = 0; l < n; l++)
            Winv[k*n + l] = sqrt(1.0 / (alpha[k] * W[l] + beta[k]));
	/*printf("D(1,    1 )     = %8f\n", Winv[0]);*/
	/*printf("D(500, 500)     = %8f\n", Winv[499]);*/
    /* sqrt(Winv) * ZtY */
    for (k = 0; k < y_inc; k++)
        for (l = 0; l < n; l++)
            y_cur[k*n+l] *= Winv[k*n+l];

	/*printf("K(1, %2d)        = %8f\n", j+1, y_cur[0]);*/
	/*printf("K(500, %2d)      = %8f\n", j+1, y_cur[499]);*/

	/*printf("xl value: %f\n", xl[0]);*/
      /* sqrt(Winv) * Z' * XL */
	for (ll = 0; ll < y_inc; ll++)
	  for ( k = 0; k < wXL; k++ )
		  for ( l = 0; l < n; l++ )
			  xl[ ll * wXL * n + k * n + l ] *= Winv[l];
	/*printf("xl value: %f\n", xl[0]);*/
	  /*printf("W(1, 1)         = %8f\n", x_cur[0]);*/
	  /*printf("W(500, 3000)    = %8f\n", x_cur[x_b*n*p-1]);*/
      
      /* 7) b = sqrt(Winv) * ZtXL * y */
	for (ll = 0; ll < y_inc; ll++)
      dgemv_("T", &n, &wXL, &ONE, &xl[ll * wXL * n], &n, &y_cur[ll * n], &iONE, &ZERO, &xl_b[ll * wXL], &iONE);
	/*printf("xl_b value: %f\n", xl_b[0]);*/

	for (ll = 0; ll < y_inc; ll++)
	{
		/*dsyrk_("L", "T", &p, &n, &ONE, &x_cur[i*p*n], &n, &ZERO, xtSx, &p);*/
		dsyrk_("L", "T", // LOWER, NO_TRANS, 
				&wXL, &n, // n, k
				&ONE, &xl[ll * wXL * n], &n, // KL KL' 
				&ZERO, &xltxl[ll * wXL * wXL], &wXL); // TL of xtSx
	}

    END_TIMING(cf->time->compute_time);
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
	    struct aiocb *aiocb_x_next_p = (struct aiocb *)aiocb_x_next_l[0];
		fgls_aio_read( aiocb_x_next_p,
				       fileno( loops_t->XR_fp ), x_next,
					   (ib + x_b) >= m ? MIN( x_b, m ) * wXR * n * sizeof(double) : MIN( x_b * wXR * n, (m - (ib + x_b)) * wXR * n ) * sizeof(double),
					   (ib + x_b) >= m ? 0 : (ib + x_b) * wXR * n * sizeof(double) );

		/*bzero( (char *)aiocb_x_next_p,  sizeof(struct aiocb) );*/
		/*aiocb_x_next_p->aio_buf = x_next;*/
		/*aiocb_x_next_p->aio_fildes = fileno( loops_t->XR_fp );*/
		/*aiocb_x_next_p->aio_nbytes = (ib + x_b) >= m ? MIN( x_b, m ) * wXR * n * sizeof(double) : MIN( x_b * wXR * n, (m - (ib + x_b)) * wXR * n ) * sizeof(double);*/
		/*aiocb_x_next_p->aio_offset = (ib + x_b) >= m ? 0 : (ib + x_b) * wXR * n * sizeof(double);*/
		/**//*aiocb_x.aio_flag = AIO_RAW;*/
		/*aio_read( aiocb_x_next_p );*/
#if VAMPIR
		VT_USER_END("READ_X");
#endif

#if VAMPIR
		VT_USER_START("WAIT_X");
#endif
		BEGIN_TIMING();
		/*aio_suspend( aiocb_x_cur_l, 1, NULL );*/
		if ( aio_suspend( aiocb_x_cur_l, 1, NULL ) != 0 )
			perror("Suspend X error\n");
	    END_TIMING(cf->time->comp_mutex_wait_time);
#if VAMPIR
		VT_USER_END("WAIT_X");
#endif
#if VAMPIR
      VT_USER_START("COMP_LOOP_X_CODE");
#endif
      BEGIN_TIMING();
      x_inc = MIN(x_b, m - ib);
	  /*mpb = p * x_b;*/
	  /*mpb_real = p * x_inc;*/

	  for ( j = 0; j < y_inc; j++ )
	  {
		  /*printf("Iter j: %d\n", j);*/
		  /* sqrt(Winv) * Z' * XR */
		  for ( k = 0; k < x_inc*wXR; k++ )
			  for ( l = 0; l < n; l++ )
				  x_copy[ k * n + l ] = x_cur[ k * n + l ] * Winv[l];
		  /*printf("W(1, 1)         = %8f\n", x_cur[0]);*/
		  /*printf("W(%d, %d)    = %8f\n", n, x_b*wXR, x_cur[x_b*n*wXR-1]);*/
		  
		  /* 7) b = WinvZtX * y */
		  /*dgemv_("T", &n, &mpb_real, &ONE, x_cur, &n, y_cur, &iONE, &ZERO, b_cur, &iONE);*/
		  /*printf("RHS(1, %2d)      = %8f\n", j+1, b_cur[0]);*/
		  /*printf("RHS(3000, %2d)   = %8f\n", j+1, b_cur[2999]);*/
		  
		  for (i = 0; i < x_inc; i++)
		  {
			  Bij = &b_cur[ j*p*x_inc   + i*p];
			  Vij = &v_cur[ j*p*p*x_inc + i*p*p];
			  /*printf("Iter i: %d (%d, %d, %d)\n", i, jb, ib, j);*/
			  /*int k;*/
			  /*printf("value: %f\n", b_cur[j * x_inc * p + i*p]);*/
			  memcpy(Bij, &xl_b[j * wXL], wXL * sizeof(double));
			  /*printf("value: %f\n", b_cur[j * x_inc * p + i*p]);*/
			  /*printf("DGEMV( T, %d, %d, %2f, %p, %d, %p, %d, %2f, %p, %d);\n",*/
			  /*n, wXR, ONE, &x_cur[i * wXR * n], n, y_cur, iONE, ZERO, b_cur[i*p + wXL], iONE);*/
			  dgemv_("T", 
					  &n, &wXR, 
					  &ONE, &x_copy[i * wXR * n], &n, &y_cur[j* n], &iONE, 
					  &ZERO, &Bij[wXL], &iONE);
			  /*printf("value: %f\n", b_cur[j * x_inc * p + i*p]);*/
			  /*printf("RHS(%2d, %2d)    = %8f\n", (ib+i)*p+1,   j+1, b_cur[(ib+i)*p]);*/
			  /*printf("RHS(%2d, %2d)    = %8f\n", (ib+i+1)*p, j+1, b_cur[(ib+i+1)*p-1]);*/
			  /*printf("Thread: %d - Iterating X(%d) - Solving\n", id, i);*/
				/* 5) W = XtZWinv * K^T */
				for( k = 0; k < wXL; k++ )
					dcopy_(&wXL, &xltxl[j*wXL*wXL + k*wXL], &iONE, &Vij[k*p], &iONE); // TL
				dgemm_("T", "N", // NO_TRANS, TRANS, 
						&wXR, &wXL, &n, // m, n, k
						&ONE, &x_copy[i * wXR * n], &n, &xl[j*wXL*n], &n, // KR KL'
						&ZERO, &Vij[wXL], &p); // xtSx BL
				dsyrk_("L", "T", 
						&wXR, &n, 
						&ONE, &x_copy[i * wXR * n], &n, 
						&ZERO, &Vij[wXL * p + wXL], &p);

				/*printf("%8f %8f %8f %8f\n", xtSx[0], xtSx[4], xtSx[8], xtSx[12]);*/
				/*printf("%8f %8f %8f %8f\n", xtSx[1], xtSx[5], xtSx[9], xtSx[13]);*/
				/*printf("%8f %8f %8f %8f\n", xtSx[2], xtSx[6], xtSx[10], xtSx[14]);*/
				/*printf("%8f %8f %8f %8f\n", xtSx[3], xtSx[7], xtSx[11], xtSx[15]);*/
				/* 8) W^-1 * y */
				/*printf("value pre dposv: %f\n", b_cur[j * x_inc * p + i*p]);*/
				dpotrf_(LOWER, &p, Vij, &p, &info);
				if (info != 0)
				{
					char err[STR_BUFFER_SIZE];
					snprintf(err, STR_BUFFER_SIZE, "dpotrf failed (info: %d)", info);
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

				/* dposv_("L", &p, &iONE, xtSx, &p, Bij, &p, &info);
				if (info != 0)
				{
				  fprintf(stderr, "Error executing dposv (y: %d, x: %d, th: %d): %d\n", jb+j, ib+i, id, info);
				  exit(-1);
				} */
				/*printf("value post dposv: %f\n", b_cur[j * x_inc * p + i*p]);*/
				/*if ( i== 0) printf("res[%d, %2d]: %f\n", ib+i, j, *b_cur);*/
		  }
	  }
	  /*printf("B(1, %d)         = %8f\n", j+1, b_cur[0]);*/
	  /*printf("B(3000, %d)      = %8f\n", j+1, b_cur[2999]);*/
      END_TIMING(cf->time->compute_time);
#if VAMPIR
      VT_USER_END("COMP_LOOP_X_CODE");
#endif
		
#if VAMPIR
      VT_USER_START("WAIT_B");
#endif
	  if ( iter > 0)
	  {
		  /*if ( aio_suspend( aiocb_b_prev_l, prev_y_inc, NULL ) != 0 ) // FIX*/
		  /*if ( aio_suspend( aiocb_b_prev_l, y_b, NULL ) != 0 ) // FIX*/
		  /*perror("Error waiting for b");*/
		  for ( k = 0; k < y_b; k++ )
		  {
			if ( aio_suspend( &aiocb_b_prev_l[k], 1, NULL ) != 0 )
			  perror("Suspend B error");
			if ( aio_suspend( &aiocb_v_prev_l[k], 1, NULL ) != 0 )
			  perror("Suspend V error");
			/*printf("%d - %d - %p - %d\n", id, aio_error((const struct aiocb*)aiocb_b_prev_l[k]), aiocb_b_prev_l[k], aiocb_b_prev_l[k]->aio_offset);*/
		  }
	  }
#if VAMPIR
      VT_USER_END("WAIT_B");
#endif

#if VAMPIR
      VT_USER_START("WRITE_B");
#endif
      bzero( (char *)aiocb_b_cur, y_inc * sizeof(struct aiocb) );
	  /*bzero( (char *)aiocb_b_cur_l[0], y_inc * sizeof(struct aiocb) );*/
	  struct aiocb *aiocb_b_cur_p;
      bzero( (char *)aiocb_v_cur, y_inc * sizeof(struct aiocb) );
	  struct aiocb *aiocb_v_cur_p;
	  for ( k = 0; k < y_inc; k++ )
	  {
		  /*aiocb_b_cur_p = &aiocb_b_cur_l[0][k];*/
		  aiocb_b_cur_p = aiocb_b_cur_l[k];
		  fgls_aio_write( aiocb_b_cur_p,
				          fileno( loops_t->B_fp ), &b_cur[ k * x_inc * p],
						  x_inc * p * sizeof(double),
						  ((jb+k) * m * p + ib * p) * sizeof(double) );
		  /*aiocb_b_cur_p->aio_reqprio = 0;*/
		  /*aiocb_b_cur_p->aio_buf = &b_cur[ k * x_inc * p];*/
		  /*aiocb_b_cur_p->aio_fildes = fileno( loops_t->B_fp );*/
		  /*aiocb_b_cur_p->aio_nbytes = x_inc * p * sizeof(double);*/
		  /*aiocb_b_cur_p->aio_offset = ((jb+k) * m * p + ib * p) * sizeof(double);*/
		  /*aiocb_b_cur_p->lio_opcode = LIO_WRITE;*/
		  /*printf("Nbytes: %d\n", aiocb_b_cur_p->aio_nbytes);*/
		  /*printf("%d - Offset: %d\n", id, aiocb_b_cur_p->aio_offset);*/
		  /*printf("%d - %p\n", id, aiocb_b_cur_l[k]);*/
		  /*aiocb_y.aio_flag = AIO_RAW;*/
		  /*if ( aio_write( aiocb_b_cur_p ) != 0 )*/
		  /*perror("Error writing b");*/

		  aiocb_v_cur_p = aiocb_v_cur_l[k];
		  fgls_aio_write( aiocb_v_cur_p,
				          fileno( loops_t->V_fp ), &v_cur[ k * x_inc * p * p],
						  x_inc * p * p * sizeof(double),
						  ((jb+k) * m * p * p + ib * p * p) * sizeof(double) );
		  /*aiocb_v_cur_p->aio_reqprio = 0;*/
		  /*aiocb_v_cur_p->aio_buf = &v_cur[ k * x_inc * p * p];*/
		  /*aiocb_v_cur_p->aio_fildes = fileno( loops_t->V_fp );*/
		  /*aiocb_v_cur_p->aio_nbytes = x_inc * p * p * sizeof(double);*/
		  /*aiocb_v_cur_p->aio_offset = ((jb+k) * m * p * p + ib * p * p) * sizeof(double);*/
		  /*if ( aio_write( aiocb_v_cur_p ) != 0 )*/
		  /*perror("Error writing v");*/
	  }
#if VAMPIR
      VT_USER_END("WRITE_B");
#endif

	  swap_aiocb( &aiocb_x_cur_l, &aiocb_x_next_l );
	  swap_buffers( &x_cur, &x_next);
	  swap_aiocb( &aiocb_b_cur_l, &aiocb_b_prev_l );
	  swap_buffers( &b_cur, &b_prev);
	  swap_aiocb( &aiocb_v_cur_l, &aiocb_v_prev_l );
	  swap_buffers( &v_cur, &v_prev);
	  iter++;
    }
    gettimeofday(&t1, NULL);
	if (id == 0)
	{
		printf("Iter (%d - %d) time: %ld ms\n", jb, jb+y_inc, get_diff_ms(&t0, &t1));
		fflush(stdout);
	}
	swap_aiocb( &aiocb_y_cur_l, &aiocb_y_next_l );
	swap_buffers( &y_cur, &y_next);
	/*prev_y_inc = y_inc;*/
  }
#if VAMPIR
      VT_USER_START("WAIT_X");
#endif
  aio_suspend( aiocb_x_cur_l, 1, NULL );
#if VAMPIR
      VT_USER_END("WAIT_X");
#endif
#if VAMPIR
      VT_USER_START("WAIT_Y");
#endif
  aio_suspend( aiocb_y_cur_l, 1, NULL );
#if VAMPIR
      VT_USER_END("WAIT_Y");
#endif
  /*aio_suspend( aiocb_b_prev_l, 1, NULL );*/
  /*if ( aio_suspend( aiocb_b_prev_l, y_inc, NULL ) != 0 )*/
#if VAMPIR
      VT_USER_START("WAIT_B");
#endif
  for ( i = 0; i < y_inc; i++ )
  {
    if ( aio_suspend( &aiocb_b_prev_l[i], 1, NULL ) != 0 )
      perror("Suspend B error");
    if ( aio_suspend( &aiocb_v_prev_l[i], 1, NULL ) != 0 )
      perror("Suspend V error");
	/*printf("END: %d - %d - %p - %d\n", id, aio_error((const struct aiocb*)aiocb_b_prev_l[i]), aiocb_b_prev_l[i], aiocb_b_prev_l[i]->aio_offset);*/
  }
#if VAMPIR
      VT_USER_END("WAIT_B");
#endif
  pthread_exit(NULL);
}

int preloop(FGLS_config_t *cf, double *Phi, double *Z, double *W) 
{
  DEF_TIMING();

  /* Threads for ooc gemms */
  pthread_t io_thread;
  pthread_t compute_thread;
  ooc_gemm_t gemm_t;
  int iret;

  /* Buffer sizes for ooc gemms */
  long int chunk_size = 1L << 28; // 256 MElems
  chunk_size = chunk_size - chunk_size % (cf->n * sizeof(double));
  int num_cols = chunk_size / (cf->n * sizeof(double));

  /* Z W Z' = Phi */
  printf("\nEigendecomposition of Phi...");
  fflush(stdout);
  
  BEGIN_TIMING();
  eigenDec( cf->n, Phi, Z, W );
  END_TIMING(cf->time->compute_time);
  printf(" Done\n");

  /* OOC gemms */
  gemm_t.in[0]  = ( double * ) malloc ( chunk_size * sizeof(double) );
  gemm_t.in[1]  = ( double * ) malloc ( chunk_size * sizeof(double) );
  gemm_t.out[0] = ( double * ) malloc ( chunk_size * sizeof(double) );
  gemm_t.out[1] = ( double * ) malloc ( chunk_size * sizeof(double) );
  if ( gemm_t.out[1] == NULL )
  {
	  fprintf(stderr, "Not enough memory for OOC gemm's\n");
	  exit(EXIT_FAILURE);
  }

  printf("Computing Z' XL...");
  fflush(stdout);

  gemm_t.Z = Z;

  gemm_t.m = cf->n;
  gemm_t.n = cf->wXL;
  gemm_t.k = cf->n;
  gemm_t.n_cols_per_buff = num_cols;
  gemm_t.fp_in  = fopen( cf->XL_path, "rb" );
  gemm_t.fp_out = fopen( cf->ZtXL_path, "wb" );
  /*fd_in  = open( cf->X_path, O_RDONLY | O_SYNC );*/
  /*fd_out = open( cf->ZtX_path, O_WRONLY | O_CREAT | O_SYNC );*/
  /*gemm_t.fp_in  = fdopen( fd_in,  "rb" );*/
  /*gemm_t.fp_out = fdopen( fd_out, "wb" );*/

  /*sem_init( &gemm_t.sem_io,   0, 0 );*/
  /*sem_init( &gemm_t.sem_comp, 0, 0 );*/

  /*iret = pthread_create(&io_thread, NULL, ooc_gemm_io, (void*)&gemm_t);*/
  /*if (iret)*/
  /*{*/
  /*fprintf(stderr, __FILE__ ": Error creating IO thread (1): %d\n", iret);*/
  /*exit(EXIT_FAILURE);*/
  /*}*/
  iret = pthread_create(&compute_thread, NULL, ooc_gemm_comp, (void*)&gemm_t);
  if (iret)
  {
    fprintf(stderr, __FILE__ ": Error creating Computation thread (1): %d\n", iret);
    exit(EXIT_FAILURE);
  }

  /*pthread_join(io_thread, NULL);*/
  pthread_join(compute_thread, NULL);

  /*printf("Closing files\n"); */
  /*printf("Closing input\n"); */
  fclose( gemm_t.fp_in );
  /*printf("Closing output\n"); */
  fclose( gemm_t.fp_out );
  /*printf("Closed\n"); */

  printf(" Done\n");

  printf("Computing Z' XR...");
  fflush(stdout);

  /*gemm_t.Z = Z;*/

  gemm_t.m = cf->n;
  gemm_t.n = cf->m * cf->wXR;
  gemm_t.k = cf->n;
  gemm_t.n_cols_per_buff = num_cols;
  gemm_t.fp_in  = fopen( cf->XR_path, "rb" );
  gemm_t.fp_out = fopen( cf->ZtXR_path, "wb" );

  /*sem_init( &gemm_t.sem_io,   0, 0 );*/
  /*sem_init( &gemm_t.sem_comp, 0, 0 );*/

  /*iret = pthread_create(&io_thread, NULL, ooc_gemm_io, (void*)&gemm_t);*/
  /*if (iret)*/
  /*{*/
  /*fprintf(stderr, __FILE__ ": Error creating IO thread (2): %d\n", iret);*/
  /*exit(EXIT_FAILURE);*/
  /*}*/
  iret = pthread_create(&compute_thread, NULL, ooc_gemm_comp, (void*)&gemm_t);
  if (iret)
  {
    fprintf(stderr, __FILE__ ": Error creating Computation thread (2): %d\n", iret);
    exit(EXIT_FAILURE);
  }

  /*pthread_join(io_thread, NULL);*/
  pthread_join(compute_thread, NULL);

  fclose( gemm_t.fp_in );
  fclose( gemm_t.fp_out );

  printf(" Done\n");

  printf("Computing Z' Y...");
  fflush(stdout);

  gemm_t.m = cf->n;
  gemm_t.n = cf->t;
  gemm_t.k = cf->n;
  gemm_t.n_cols_per_buff = num_cols;
  gemm_t.fp_in  = fopen( cf->Y_path, "rb" );
  gemm_t.fp_out = fopen( cf->ZtY_path, "wb" );

  /*sem_init( &gemm_t.sem_io,   0, 0 );*/
  /*sem_init( &gemm_t.sem_comp, 0, 0 );*/

  /*iret = pthread_create(&io_thread, NULL, ooc_gemm_io, (void*)&gemm_t);*/
  /*if (iret)*/
  /*{*/
  /*fprintf(stderr, __FILE__ ": Error creating IO thread (3): %d\n", iret);*/
  /*exit(EXIT_FAILURE);*/
  /*}*/
  iret = pthread_create(&compute_thread, NULL, ooc_gemm_comp, (void*)&gemm_t);
  if (iret)
  {
    fprintf(stderr, __FILE__ ": Error creating Computation thread (3): %d\n", iret);
    exit(EXIT_FAILURE);
  }

  /*pthread_join(io_thread, NULL);*/
  pthread_join(compute_thread, NULL);

  fclose( gemm_t.fp_in );
  fclose( gemm_t.fp_out );

  printf(" Done\n");

  free(gemm_t.in[0]);
  free(gemm_t.in[1]);
  free(gemm_t.out[0]);
  free(gemm_t.out[1]);

  sem_destroy( &gemm_t.sem_io );
  sem_destroy( &gemm_t.sem_comp );

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

	work   = (double *) malloc ( lwork * sizeof(double) );
	iwork  = (int *)    malloc ( liwork * sizeof(int) );
	isuppz = (int *)    malloc ( 2 * n * sizeof(int) );

	if ( isuppz == NULL )
	{
		fprintf( stderr, __FILE__ "Not enough memory\n" );
		exit( EXIT_FAILURE );
	}

	dsyevr_("V", "A", "L", &n, Phi, &n, 
			&ddummy, &ddummy, &idummy, &idummy, &ddummy, 
			&nCompPairs, W, Z, &n, isuppz, 
			work, &lwork, iwork, &liwork, &info);
	if (info != 0)
	{
		char err[STR_BUFFER_SIZE];
		snprintf(err, STR_BUFFER_SIZE, "dsyevr failed (info: %d)", info);
		error_msg(err, 1);
	}
}

/*void* ooc_gemm_io( void *in ) */
/*{*/
/*DEF_TIMING();*/

/*ooc_gemm_t *gemm_t = ( ooc_gemm_t* )in;*/
/*#if TIMING*/
/*FGLS_eigen_t *cf = &FGLS_eigen_config;*/
/*#endif*/
/*double *in_cur   = gemm_t->in[0];*/
/*double *in_next  = gemm_t->in[1];*/
/*double *out_prev = gemm_t->out[0];*/
/*double *out_cur  = gemm_t->out[1];*/

/*int m = gemm_t->m,*/
/*n = gemm_t->n,*/
/*k = gemm_t->k,*/
/*cols_per_buff = gemm_t->n_cols_per_buff;*/
/*long int max_elems = cols_per_buff * m;*/

/*int i;*/

/**//*printf("io m: %d\n", m);*/
/**//*printf("io n: %d\n", n);*/
/**//*printf("io k: %d\n", k);*/
/*BEGIN_TIMING();*/
/*sync_read( in_cur, gemm_t->fp_in, MIN( max_elems, k * n ), 0 );*/
/*END_TIMING(cf->time->io_time);*/

/*printf("buf[1]: %f\n", in_cur[1]);*/

/*sem_post( &gemm_t->sem_comp );*/

/**//*printf("m: %d\n", m);*/
/**//*printf("n: %d\n", n);*/
/**//*printf("k: %d\n", k);*/
/**//*printf("n_cols: %d\n", cols_per_buff);*/
/**//*printf("chunk:: %d\n", chunk_size);*/
/*for ( i = 0; i < n; i += cols_per_buff ) */
/*{*/
/*swap_buffers(&in_cur, &in_next);*/

/*BEGIN_TIMING();*/
/*sync_read( in_cur, gemm_t->fp_in, */
/*i + cols_per_buff > n ? 0 : MIN( max_elems, ( n - ( i + cols_per_buff ) ) * k ), */
/*(i + cols_per_buff) * k );*/
/*END_TIMING(cf->time->io_time);*/

/*sem_post( &gemm_t->sem_comp );*/

/*BEGIN_TIMING();*/
/*sem_wait( &gemm_t->sem_io );*/
/*END_TIMING(cf->time->io_mutex_wait_time);*/

/*BEGIN_TIMING();*/
/**//*printf("Writing: %f\n", out_prev[0]);*/
/**//*printf("Writing: %f\n", out_prev[MIN( max_elems, (n - i) * m )-1]);*/
/*sync_write( out_prev, gemm_t->fp_out, MIN( max_elems, (n - i) * m ), i * m);*/
/**//*fflush( gemm_t->fp_out );*/
/*END_TIMING(cf->time->io_time);*/
/*swap_buffers( &out_prev, &out_cur );*/
/*}*/
/*pthread_exit(NULL);*/
/*}*/

void* ooc_gemm_comp( void *in ) 
{
  ooc_gemm_t *gemm_t = ( ooc_gemm_t *)in;

  double *in_cur   = gemm_t->in[0];
  double *in_next  = gemm_t->in[1];
  double *out_prev = gemm_t->out[0];
  double *out_cur  = gemm_t->out[1];

  int m = gemm_t->m,
	  n = gemm_t->n,
	  k = gemm_t->k,
	  n_cols_per_buff = gemm_t->n_cols_per_buff;
  long int max_elems = n_cols_per_buff * m;

  double ONE  = 1.0;
  double ZERO = 0.0;

  struct aiocb aiocb_in_cur,  aiocb_in_next,
			   aiocb_out_prev, aiocb_out_cur;

  struct aiocb **aiocb_in_cur_l,
			   **aiocb_in_next_l,
			   **aiocb_out_prev_l,
			   **aiocb_out_cur_l;

  aiocb_in_cur_l   = (struct aiocb **) malloc (sizeof(struct aiocb *));
  aiocb_in_next_l  = (struct aiocb **) malloc (sizeof(struct aiocb *));
  aiocb_out_prev_l = (struct aiocb **) malloc (sizeof(struct aiocb *));
  aiocb_out_cur_l  = (struct aiocb **) malloc (sizeof(struct aiocb *));

  aiocb_in_cur_l[0]   = &aiocb_in_cur;
  aiocb_in_next_l[0]  = &aiocb_in_next;
  aiocb_out_prev_l[0] = &aiocb_out_prev;
  aiocb_out_cur_l[0]  = &aiocb_out_cur;

  /*bzero( (char *)&aiocb_in_cur,   sizeof(struct aiocb) );*/
  /*bzero( (char *)&aiocb_in_next,  sizeof(struct aiocb) );*/

  fgls_aio_read( &aiocb_in_cur,
		         fileno( gemm_t->fp_in ), in_cur,
				 MIN( max_elems, k * n ) * sizeof(double), 0);
  /*aiocb_in_cur.aio_fildes = fileno( gemm_t->fp_in );*/
  /*aiocb_in_cur.aio_buf = in_cur;*/
  /*aiocb_in_cur.aio_nbytes = MIN( max_elems, k * n ) * sizeof(double);*/
  /*aiocb_in_cur.aio_offset = 0;*/
  /*aio_read( &aiocb_in_cur );*/

  int cur_n;
  int i;
  for ( i = 0; i < n; i += n_cols_per_buff ) 
  {
	struct aiocb *aiocb_in_next_p = (struct aiocb *)aiocb_in_next_l[0];
	fgls_aio_read( aiocb_in_next_p,
			       fileno( gemm_t->fp_in ), in_next,
				   i + n_cols_per_buff > n ? 0 : MIN( max_elems, ( n - ( i + n_cols_per_buff ) ) * k ) * sizeof(double),
				   (i + n_cols_per_buff) * k * sizeof(double) );

	/*bzero( (char *)aiocb_in_next_p,  sizeof(struct aiocb) );*/
	/*aiocb_in_next_p->aio_fildes = fileno( gemm_t->fp_in );*/
	/*aiocb_in_next_p->aio_buf = in_next;*/
	/*aiocb_in_next_p->aio_nbytes = i + n_cols_per_buff > n ? 0 : MIN( max_elems, ( n - ( i + n_cols_per_buff ) ) * k ) * sizeof(double);*/
	/*aiocb_in_next_p->aio_offset = (i + n_cols_per_buff) * k * sizeof(double);*/
	/*aio_read( aiocb_in_next_p );*/

	if ( aio_suspend( aiocb_in_cur_l, 1, NULL ) != 0 )
		perror("Suspend In error\n");

	cur_n = MIN( n_cols_per_buff, (n - i) );
    dgemm_("T", "N", &m, &cur_n, &k, &ONE, gemm_t->Z, &m, in_cur, &m, &ZERO, out_cur, &m);

	if ( i > 0)
	{
		if ( aio_suspend( aiocb_out_prev_l, 1, NULL ) != 0 )
		  perror("Suspend Out error");
	}

	struct aiocb *aiocb_out_cur_p = (struct aiocb *)aiocb_out_cur_l[0];
	fgls_aio_write( aiocb_out_cur_p,
			        fileno( gemm_t->fp_out ), out_cur,
					MIN( max_elems, (n - i) * m ) * sizeof(double),
					i * m * sizeof(double) );
	/*bzero( (char *)aiocb_out_cur_p, sizeof(struct aiocb) );*/

	/**//*aiocb_out_cur_p = aiocb_out_cur_l;*/
	/*aiocb_out_cur_p->aio_buf = out_cur;*/
	/*aiocb_out_cur_p->aio_fildes = fileno( gemm_t->fp_out );*/
	/*aiocb_out_cur_p->aio_nbytes = MIN( max_elems, (n - i) * m ) * sizeof(double);*/
	/*aiocb_out_cur_p->aio_offset = i * m * sizeof(double);*/

	/*if ( aio_write( aiocb_out_cur_p ) != 0 )*/
	/*perror("Error writing Out");*/

	swap_aiocb( &aiocb_in_cur_l, &aiocb_in_next_l );
	swap_aiocb( &aiocb_out_prev_l, &aiocb_out_cur_l );
    swap_buffers( &in_cur,  &in_next );
    swap_buffers( &out_prev, &out_cur );
  }

  if ( aio_suspend( aiocb_in_cur_l, 1, NULL ) != 0 )
	perror("Suspend In error\n");
  if ( aio_suspend( aiocb_out_prev_l, 1, NULL ) != 0 )
    perror("Suspend Out error");

  pthread_exit(NULL);
}
