#include <stdio.h>
#include <stdlib.h>

#include <sys/types.h>
#include <sys/stat.h>
/*#include <fcntl.h>*/

#include <pthread.h>
#include <semaphore.h>

#include "options.h"
#include "blas.h"
#include "lapack.h"

#include "timing.h"
#include "io.h"
#include "fgls_eigen.h"

typedef struct 
{
	double *Z;
	FILE *fp_in;
	FILE *fp_out;
	int m, n, k;
	long int n_cols_per_buff;
	double *in[2];
	double *out[2];
	sem_t sem_io;
	sem_t sem_comp;
} ooc_gemm_t;


/* 
 * Eigendecomposition of Phi
 *
 * Z W Z' = Phi
 */
void eigenDec(int n, double *Phi, double *Z, double *W)
{
	int nb = 192;
	int idummy, nCompPairs, *isuppz, *iwork, info,
		lwork = n * (nb + 6),
		liwork = 10 * n;
	double ddummy = -1.0, *work;

	work   = (double *) malloc ( lwork * sizeof(double) );
	iwork  = (int *)    malloc ( liwork * sizeof(int) );
	isuppz = (int *)    malloc ( 2 * n * sizeof(int) );

	dsyevr_("V", "A", "L", &n, Phi, &n, 
			&ddummy, &ddummy, &idummy, &idummy, &ddummy, 
			&nCompPairs, W, Z, &n, isuppz, 
			work, &lwork, iwork, &liwork, &info);

	if (info != 0) 
	{
		fprintf(stderr, __FILE__ ": Error while decomposing Phi\n");
		exit(-1);
	}
}

void* ooc_gemm_io( void *in ) 
{
  DEF_TIMING();

  ooc_gemm_t *gemm_t = ( ooc_gemm_t* )in;
  FGLS_eigen_t *cf = &FGLS_eigen_config;

  double *in_cur   = gemm_t->in[0];
  double *in_next  = gemm_t->in[1];
  double *out_prev = gemm_t->out[0];
  double *out_cur  = gemm_t->out[1];

  int m = gemm_t->m,
	  n = gemm_t->n,
	  k = gemm_t->k,
	  cols_per_buff = gemm_t->n_cols_per_buff;
  long int max_elems = cols_per_buff * m;

  int i;

  /*printf("io m: %d\n", m);*/
  /*printf("io n: %d\n", n);*/
  /*printf("io k: %d\n", k);*/
  BEGIN_TIMING();
  sync_read( in_cur, gemm_t->fp_in, MIN( max_elems, k * n ), 0 );
  END_TIMING(cf->time->io_time);

  sem_post( &gemm_t->sem_comp );

  /*printf("m: %d\n", m);*/
  /*printf("n: %d\n", n);*/
  /*printf("k: %d\n", k);*/
  /*printf("n_cols: %d\n", cols_per_buff);*/
  /*printf("chunk:: %d\n", chunk_size);*/
  for ( i = 0; i < n; i += cols_per_buff ) 
  {
    swap_buffers(&in_cur, &in_next);
    
    BEGIN_TIMING();
    sync_read( in_cur, gemm_t->fp_in, 
			i + cols_per_buff > n ? 0 : MIN( max_elems, ( n - ( i + cols_per_buff ) ) * k ), 
			(i + cols_per_buff) * k );
    END_TIMING(cf->time->io_time);

    sem_post( &gemm_t->sem_comp );
    
    BEGIN_TIMING();
    sem_wait( &gemm_t->sem_io );
    END_TIMING(cf->time->io_mutex_wait_time);

    BEGIN_TIMING();
	/*printf("Writing: %f\n", out_prev[0]);*/
	/*printf("Writing: %f\n", out_prev[MIN( max_elems, (n - i) * m )-1]);*/
    sync_write( out_prev, gemm_t->fp_out, MIN( max_elems, (n - i) * m ), i * m);
	/*fflush( gemm_t->fp_out );*/
    END_TIMING(cf->time->io_time);
    swap_buffers( &out_prev, &out_cur );
  }
  pthread_exit(NULL);
}

void* ooc_gemm_comp( void *in ) 
{
  DEF_TIMING();

  ooc_gemm_t *gemm_t = ( ooc_gemm_t *)in;
  FGLS_eigen_t *cf = &FGLS_eigen_config;

  double *in_cur   = gemm_t->in[0];
  double *in_next  = gemm_t->in[1];
  double *out_cur  = gemm_t->out[0];
  double *out_next = gemm_t->out[1];

  int m = gemm_t->m,
	  n = gemm_t->n,
	  k = gemm_t->k,
	  n_cols_per_buff = gemm_t->n_cols_per_buff;

  double ONE  = 1.0;
  double ZERO = 0.0;

  int cur_n;
  int i;
  /*printf("m: %d\n", m);*/
  /*printf("n: %d\n", n);*/
  /*printf("k: %d\n", k);*/
  for ( i = 0; i < n; i += n_cols_per_buff ) 
  {
	  /*printf("computing from col: %d\n", i);*/
    BEGIN_TIMING();
    sem_wait( &gemm_t->sem_comp );
    END_TIMING(cf->time->comp_mutex_wait_time);
	/*printf("computing from col(2): %d\n", i);*/

    BEGIN_TIMING();
	cur_n = MIN( n_cols_per_buff, (n - i) );
    dgemm_("T", "N", &m, &cur_n, &k, &ONE, gemm_t->Z, &m, in_cur, &m, &ZERO, out_cur, &m);
    END_TIMING(cf->time->compute_time);

    sem_post( &gemm_t->sem_io );
    swap_buffers( &in_cur,  &in_next );
    swap_buffers( &out_cur, &out_next );
  }
  pthread_exit(NULL);
}

int preloop(double *Phi, double *Z, double *W) 
{
  DEF_TIMING();

  int iret;
  void *vret;
  pthread_t io_thread;
  pthread_t compute_thread;
  ooc_gemm_t gemm_t;
  long GB = 1L << 24;

  FGLS_eigen_t *cf = &FGLS_eigen_config;

  long int chunk_size = GB - GB % (cf->n * sizeof(double));
  int num_cols = chunk_size / (cf->n * sizeof(double));
  /*printf("Chunk: %ld MB\n", chunk_size >> 20);*/
  /*printf("numcols: %ld\n", num_cols);*/

  int fd_in, fd_out;

  /*printf("chunk size: %d\n", chunk_size);*/
  /*printf("Num cols: %d\n", num_cols);*/
  /*printf("n: %d\n", cf->n);;*/

  /*printf("Z[0]: %f\n", Z[0]);*/

  /* Z W Z' = Phi */
  printf("\nEigendecomposition of Phi...");
  fflush(stdout);
  BEGIN_TIMING();
  eigenDec( cf->n, Phi, Z, W );
  END_TIMING(cf->time->compute_time);
  printf(" Done\n");

  gemm_t.in[0]  = ( double* ) malloc ( chunk_size * sizeof(double) );
  gemm_t.in[1]  = ( double* ) malloc ( chunk_size * sizeof(double) );
  gemm_t.out[0] = ( double* ) malloc ( chunk_size * sizeof(double) );
  gemm_t.out[1] = ( double* ) malloc ( chunk_size * sizeof(double) );
  if (gemm_t.out[1] == NULL)
  {
	  fprintf(stderr, "Not enough memory for OOC gemm's\n");
	  exit(EXIT_FAILURE);
  }

  /* OOC gemms */
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
  sem_init( &gemm_t.sem_io,   0, 0 );
  sem_init( &gemm_t.sem_comp, 0, 0 );

  iret = pthread_create(&io_thread, NULL, ooc_gemm_io, (void*)&gemm_t);
  if (iret)
  {
    fprintf(stderr, __FILE__ ": Error creating IO thread (1): %d\n", iret);
    exit(EXIT_FAILURE);
  }
  iret = pthread_create(&compute_thread, NULL, ooc_gemm_comp, (void*)&gemm_t);
  if (iret)
  {
    fprintf(stderr, __FILE__ ": Error creating Computation thread (1): %d\n", iret);
    exit(EXIT_FAILURE);
  }

  pthread_join(io_thread, &vret);
  pthread_join(compute_thread, &vret);

  /*printf("Closing files\n"); */
  /*printf("Closing input\n"); */
  fclose( gemm_t.fp_in );
  /*printf("Closing output\n"); */
  fclose( gemm_t.fp_out );
  /*printf("Closed\n"); */

  printf(" Done\n");

  printf("Computing Z' XR...");
  fflush(stdout);

  /*chunk_size = MIN ( 1L<<19 - (1L << 19 % cf->n), cf->m * cf->wXR * cf->n );*/
  /*num_cols = chunk_size / cf->n;*/

  gemm_t.Z = Z;
  /*gemm_t.in[0]  = ( double* ) malloc ( chunk_size * sizeof(double) );*/
  /*gemm_t.in[1]  = ( double* ) malloc ( chunk_size * sizeof(double) );*/
  /*gemm_t.out[0] = ( double* ) malloc ( chunk_size * sizeof(double) );*/
  /*gemm_t.out[1] = ( double* ) malloc ( chunk_size * sizeof(double) );*/

  gemm_t.m = cf->n;
  gemm_t.n = cf->m * cf->wXR;
  gemm_t.k = cf->n;
  gemm_t.n_cols_per_buff = num_cols;
  gemm_t.fp_in  = fopen( cf->XR_path, "rb" );
  gemm_t.fp_out = fopen( cf->ZtXR_path, "wb" );

  sem_init( &gemm_t.sem_io,   0, 0 );
  sem_init( &gemm_t.sem_comp, 0, 0 );

  iret = pthread_create(&io_thread, NULL, ooc_gemm_io, (void*)&gemm_t);
  if (iret)
  {
    fprintf(stderr, __FILE__ ": Error creating IO thread (2): %d\n", iret);
    exit(EXIT_FAILURE);
  }
  iret = pthread_create(&compute_thread, NULL, ooc_gemm_comp, (void*)&gemm_t);
  if (iret)
  {
    fprintf(stderr, __FILE__ ": Error creating Computation thread (2): %d\n", iret);
    exit(EXIT_FAILURE);
  }

  pthread_join(io_thread, &vret);
  pthread_join(compute_thread, &vret);

  /*printf("Closing files\n"); */
  /*printf("Closing input\n"); */
  fclose( gemm_t.fp_in );
  /*printf("Closing output\n"); */
  fclose( gemm_t.fp_out );
  /*printf("Closed\n"); */

  printf(" Done\n");

  /*chunk_size = MIN ( 1L<<19 - (1L << 19 % cf->n), cf->t * cf->n );*/
  /*num_cols = chunk_size / cf->n;*/

  gemm_t.m = cf->n;
  gemm_t.n = cf->t;
  gemm_t.k = cf->n;
  gemm_t.n_cols_per_buff = num_cols;
  gemm_t.fp_in  = fopen( cf->Y_path, "rb" );
  gemm_t.fp_out = fopen( cf->ZtY_path, "wb" );

  sem_init( &gemm_t.sem_io,   0, 0 );
  sem_init( &gemm_t.sem_comp, 0, 0 );

  printf("Computing Z' Y...");
  fflush(stdout);
  iret = pthread_create(&io_thread, NULL, ooc_gemm_io, (void*)&gemm_t);
  if (iret)
  {
    fprintf(stderr, __FILE__ ": Error creating IO thread (3): %d\n", iret);
    exit(EXIT_FAILURE);
  }
  iret = pthread_create(&compute_thread, NULL, ooc_gemm_comp, (void*)&gemm_t);
  if (iret)
  {
    fprintf(stderr, __FILE__ ": Error creating Computation thread (3): %d\n", iret);
    exit(EXIT_FAILURE);
  }

  pthread_join(io_thread, &vret);
  pthread_join(compute_thread, &vret);

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
