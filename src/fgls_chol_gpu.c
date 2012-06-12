#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <sys/types.h>
#include <sys/stat.h>

#include <sys/time.h>
#include <time.h>

#include <aio.h>

#include <cuda_runtime_api.h>
#include <cublas_v2.h>

#include "blas.h"
#include "lapack.h"
#include "options.h"
#include "common.h"
#include "io.h"
#include "timing.h"
#include "fgls_chol.h"

#ifdef VTRACE
    #include "vt_user.h"
#endif

void sync_gpus(int ngpus)
{
    int igpu = 0;
    for(igpu = 0 ; igpu < ngpus ; ++igpu) {
        cudaSetDevice(igpu);
        cudaStreamSynchronize(0);
    }
}

void start_section(struct timeval* t_start,
                   int ngpus,
                   const char* vt_id,
                   const char* text, ...)
{
#ifdef FGLS_GPU_SERIAL
    sync_gpus(ngpus);
#endif
#ifdef TIMING
    read_clock(t_start);
#endif
#if defined(DEBUG) || defined(TIMING)
    va_list argp;
    va_start(argp, text);
    vprintf(text, argp);
    fflush(stdout);
#endif

#ifdef VTRACE
    VT_USER_START(vt_id);
#endif
}

void end_section(struct timeval* t_start,
                 int ngpus,
                 const char* vt_id)
{
#ifdef FGLS_GPU_SERIAL
    sync_gpus(ngpus);
#endif
#ifdef VTRACE
    VT_USER_END(vt_id);
#endif
#ifdef TIMING
    struct timeval t_end;
    read_clock(&t_end);
    int dt = elapsed_time(t_start, &t_end);
    printf("done in %fs (%fms, %dus)\n", dt*0.000001, dt*0.001, dt);
    fflush(stdout);
#elif defined (DEBUG)
    printf("done\n");
    fflush(stdout);
#endif
}

#define START_SECTION(VT_ID, MSG) start_section(&t_start, ngpus, VT_ID, MSG "... ")
#define START_SECTION2(VT_ID, MSG, ...) start_section(&t_start, ngpus, VT_ID, MSG "... ", __VA_ARGS__)
#define END_SECTION(VT_ID) end_section(&t_start, ngpus, VT_ID)

/*
 * Cholesky-based solution of the Feasible Generalized Least-Squares problem
 */
int fgls_chol_gpu(int n, int p, int m, int t, int wXL, int wXR,
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

    /* Out-of-core operands (triple and doublebuffering buffers) */
    double *X[3];
    double *Y[2];
    double *B[2];
    double *V[2];
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
    int i, j, k; // size_t
    int nn = n * n; // size_t
    char numths_str[STR_BUFFER_SIZE];

    // For measuring times.
    struct timeval t_start;

#ifdef DEBUG
    /* Checking the input arguments */
    printf("n: %d\np: %d\nm: %d\nt: %d\nwXL: %d\nwXR: %d\n", n, p, m, t, wXL, wXR);
    printf("x_b: %d\ny_b: %d\nnths: %d\n", x_b, y_b, num_threads);
    printf("Phi: %s\nh2: %s\ns2: %s\n", Phi_path, h2_path, sigma2_path);
    printf("XL: %s\nXR: %s\nY: %s\n", XL_path, XR_path, Y_path);
    printf("B: %s\nV: %s\n", B_path, V_path);
#endif
#ifdef FGLS_GPU_SERIAL
    printf("Running the serialized GPU version. (Due to debug mode or FGLS_GPU_SERIAL option.)\n");
#else
    printf("Running the GPU version\n");
#endif

    if(t != 1) {
        printf("The GPU version doesn't support t != 1\n");
        exit(1);
        // Note/TODO: to make it work, you need to fix the "fetching next Xr from disk" code part to restart
        // from offset zero at some point.
    }

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

    start_section(&t_start, 0, "GPU_init", "Initializing the GPU(s)... ");

    //GPU: Initializing the GPU
    int ngpus = 0, igpu = 0;
    cublasHandle_t cu_handle;
    cudaError_t cu_error, cu_error2;
    cublasStatus_t cu_status;
    if((cu_status = cublasCreate(&cu_handle)) != CUBLAS_STATUS_SUCCESS) {
        char err[STR_BUFFER_SIZE];
        snprintf(err, STR_BUFFER_SIZE, "cublasCreate() failed (info: %d)", cu_status);
        error_msg(err, 1);
    }

    if((cu_error = cudaGetDeviceCount(&ngpus)) != CUBLAS_STATUS_SUCCESS) {
        char err[STR_BUFFER_SIZE];
        snprintf(err, STR_BUFFER_SIZE, "Can't get the cuda device count. Are there any? (info: %d)", cu_error);
        error_msg(err, 1);
    }
    END_SECTION("GPU_init");

    printf("Using %d GPUs\n", ngpus);

    // Create two streams for each GPU: one computation stream and one data transfer stream, so those two can work in parallel.
    cudaStream_t* cu_trans_streams = fgls_malloc(ngpus*sizeof(cudaStream_t));
    cudaStream_t* cu_comp_streams = fgls_malloc(ngpus*sizeof(cudaStream_t));

    for(igpu = 0 ; igpu < ngpus ; ++igpu) {
        cudaSetDevice(igpu);
        if((cu_error = cudaStreamCreate(cu_trans_streams + igpu)) != cudaSuccess) {
            char err[STR_BUFFER_SIZE];
            snprintf(err, STR_BUFFER_SIZE, "Something is wrong with your GPUs: couldn't create a transfer stream on GPU %d (info: %d)", igpu, cu_error);
            error_msg(err, 1);
        }
        if((cu_error = cudaStreamCreate(cu_comp_streams + igpu)) != cudaSuccess) {
            char err[STR_BUFFER_SIZE];
            snprintf(err, STR_BUFFER_SIZE, "Something is wrong with your GPUs: couldn't create a computation stream on GPU %d (info: %d)", igpu, cu_error);
            error_msg(err, 1);
        }
    }

    // GPU: We can already allocate GPU space for L aswell as for the streamed Xrs here.
    void** L_gpus = fgls_malloc(ngpus*sizeof(void*));
    size_t L_gpu_bytes = (size_t)cf.n * cf.n * sizeof(double);
    START_SECTION2("GPU_alloc_L", "Allocating Memory for L on each GPU: %ld bytes (%g MB)", L_gpu_bytes, L_gpu_bytes/1024.0/1024.0);
    for(igpu = 0 ; igpu < ngpus ; igpu++) {
        cudaSetDevice(igpu);
        if((cu_error = cudaMalloc(L_gpus + igpu, L_gpu_bytes)) != cudaSuccess) {
            char err[STR_BUFFER_SIZE];
            snprintf(err, STR_BUFFER_SIZE, "Not enough memory to allocate %ld bytes for L on GPU %d (info: %d)", L_gpu_bytes, igpu, cu_error);
            error_msg(err, 1);
        }
    }
    END_SECTION("GPU_alloc_L");

    double** Xr_alpha_gpus = fgls_malloc(ngpus*sizeof(double*));
    double** Xr_beta_gpus = fgls_malloc(ngpus*sizeof(double*));
    size_t Xr_gpu_bytes = (size_t)cf.x_b * cf.wXR * cf.n * sizeof(double);
    START_SECTION2("GPU_alloc_Xr", "Allocating Memory for Xr on each GPU: %ldx3 bytes (%gx3 MB)", Xr_gpu_bytes, Xr_gpu_bytes/1024.0/1024.0);
    for(igpu = 0 ; igpu < ngpus ; igpu++) {
        cudaSetDevice(igpu);
        if((cu_error = cudaMalloc((void**)(Xr_alpha_gpus + igpu), Xr_gpu_bytes)) != cudaSuccess) {
            char err[STR_BUFFER_SIZE];
            snprintf(err, STR_BUFFER_SIZE, "Not enough memory to allocate %ld bytes for Xr on GPU %d (info: %d)", Xr_gpu_bytes, igpu, cu_error);
            error_msg(err, 1);
        }
        if((cu_error = cudaMalloc((void**)(Xr_beta_gpus + igpu), Xr_gpu_bytes)) != cudaSuccess) {
            char err[STR_BUFFER_SIZE];
            snprintf(err, STR_BUFFER_SIZE, "Not enough memory to allocate %ld bytes for Xr on GPU %d (info: %d)", Xr_gpu_bytes, igpu, cu_error);
            error_msg(err, 1);
        }
    }
    END_SECTION("GPU_alloc_Xr");

    START_SECTION2("CPU_alloc", "Allocating a lot (~%g MB) of memory on the CPU", sizeof(double)*
            (cf.n*cf.n*2
            +cf.t*2
            +cf.wXL*cf.n*2
            +cf.wXL*2
            +cf.wXL*cf.wXL*2
            +cf.n
            +NUM_BUFFERS_PER_THREAD*(cf.x_b*cf.wXR*cf.n + cf.n + cf.x_b*cf.p + cf.x_b*cf.p*cf.p))/1024.0/1024.0);

    /* Memory allocation */
    // In-core
    Phi   = ( double * ) fgls_malloc ( (size_t)cf.n * cf.n * sizeof(double) );
    M     = ( double * ) fgls_malloc ( (size_t)cf.n * cf.n * sizeof(double) );
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
    // Note/TODO: allocate several smallr, non-portable blocks for each device?
    cu_error  = cudaHostAlloc((void**)(&X[0]), (size_t)cf.x_b * cf.wXR * cf.n * sizeof(double), cudaHostAllocPortable);
    cu_error2 = cudaHostAlloc((void**)(&X[1]), (size_t)cf.x_b * cf.wXR * cf.n * sizeof(double), cudaHostAllocPortable);
    if(cu_error != cudaSuccess || cu_error2 != cudaSuccess) {
        char err[STR_BUFFER_SIZE];
        snprintf(err, STR_BUFFER_SIZE, "Not enough memory to allocate 2 * %ld page-locked bytes for Xr on CPU (info: %d)", (size_t)cf.x_b * cf.wXR * cf.n * sizeof(double), cu_error);
        error_msg(err, 1);
    }

    Y[0] = ( double * ) fgls_malloc ( (size_t) cf.n * sizeof(double) );
    Y[1] = ( double * ) fgls_malloc ( (size_t) cf.n * sizeof(double) );
    B[0] = ( double * ) fgls_malloc ( (size_t) cf.x_b * cf.p * sizeof(double) );
    B[1] = ( double * ) fgls_malloc ( (size_t) cf.x_b * cf.p * sizeof(double) );
    V[0] = ( double * ) fgls_malloc ( (size_t) cf.x_b * cf.p * cf.p * sizeof(double) );
    V[1] = ( double * ) fgls_malloc ( (size_t) cf.x_b * cf.p * cf.p * sizeof(double) );
    END_SECTION("CPU_alloc");

    START_SECTION("READ_files", "Reading data for Phi, h, sigma and XL");
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
    END_SECTION("READ_files");

    /* Files and pointers for out-of-core */
    XR_fp = fopen( cf.XR_path, "rb");
    Y_fp = fopen( cf.Y_path, "rb");
    B_fp = fopen( cf.B_path, "wb");
    V_fp = fopen( cf.V_path, "wb");
    FILE* LXR_fp = fopen( "LXR.out", "wb");

    int cpu_gpu_buff = 0; // == a & b
    int hdd_cpu_buff = 1; // == c

    double *y_cur  = Y[0];
    double *y_next = Y[1];
    double *b_cur  = B[0];
    double *b_prev = B[1];
    double *v_cur  = V[0];
    double *v_prev = V[1];

    /* Asynchronous IO data structures */
    struct aiocb aiocb_x,
                 aiocb_y_cur,  aiocb_y_next,
                 aiocb_b_prev, aiocb_b_cur,
                 aiocb_v_prev, aiocb_v_cur;

    struct aiocb *aiocb_x_l[1] = {&aiocb_x};
    const struct aiocb ** aiocb_y_cur_l, // = { &aiocb_y_cur },
                       ** aiocb_y_next_l,// = { &aiocb_y_next },
                       ** aiocb_b_prev_l,// = { aiocb_b_prev },
                       ** aiocb_b_cur_l, // = { aiocb_b_cur };
                       ** aiocb_v_prev_l,// = { aiocb_v_prev },
                       ** aiocb_v_cur_l; // = { aiocb_v_cur };

    aiocb_y_cur_l  = (const struct aiocb **) fgls_malloc (sizeof(struct aiocb *));
    aiocb_y_next_l = (const struct aiocb **) fgls_malloc (sizeof(struct aiocb *));
    aiocb_b_prev_l = (const struct aiocb **) fgls_malloc (sizeof(struct aiocb *));
    aiocb_b_cur_l  = (const struct aiocb **) fgls_malloc (sizeof(struct aiocb *));
    aiocb_v_prev_l = (const struct aiocb **) fgls_malloc (sizeof(struct aiocb *));
    aiocb_v_cur_l  = (const struct aiocb **) fgls_malloc (sizeof(struct aiocb *));

    aiocb_y_cur_l[0]  = &aiocb_y_cur;
    aiocb_y_next_l[0] = &aiocb_y_next;
    aiocb_b_prev_l[0] = &aiocb_b_prev;
    aiocb_b_cur_l[0]  = &aiocb_b_cur;
    aiocb_v_prev_l[0] = &aiocb_v_prev;
    aiocb_v_cur_l[0]  = &aiocb_v_cur;

    START_SECTION("READ_Y", "Starting read of the 0th Y from disk");
    /* Read first Y */
    fgls_aio_read( &aiocb_y_cur,
                   fileno( Y_fp ), y_cur,
                   (size_t)n * sizeof(double), 0 );
    END_SECTION("READ_Y");

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

        START_SECTION2("READ_Y", "Starting read of the %dth Y from disk", j+1);
        /* Read next Y */
        struct aiocb *aiocb_y_next_p = (struct aiocb *)aiocb_y_next_l[0];
        fgls_aio_read( aiocb_y_next_p,
                       fileno( Y_fp ), y_next,
                       j + 1 >= t ? 0 : (size_t)n * sizeof(double),
                       j + 1 >= t ? 0 : (off_t)(j+1) * n * sizeof(double) );
        END_SECTION("READ_Y");

#ifdef VTRACE
        VT_USER_START("COMP_J");
#endif
        START_SECTION("COMP_M", "Computing matrix M");
        /* M := sigma * ( h^2 Phi - (1 - h^2) I ) */
        memcpy( M, Phi, (size_t)n * n * sizeof(double) );
        alpha = h[j] * sigma[j];
        beta  = (1 - h[j]) * sigma[j];
        dscal_(&nn, &alpha, M, &iONE);
        for ( i = 0; i < n; i++ )
            M[i*n + i] = M[i*n + i] + beta;
        END_SECTION("COMP_M");

        /* L * L' = M */
        START_SECTION("COMP_L", "Factorizing M into L");
        /*struct timeval t0, t1;*/
        /*gettimeofday(&t0, NULL);*/
        dpotrf_(LOWER, &n, M, &n, &info);
        END_SECTION("COMP_L");
        /*gettimeofday(&t1, NULL);*/
        /*printf("Chol:  %ld msecs\n", get_diff_ms( &t0, &t1 ));*/
        if (info != 0)
        {
            char err[STR_BUFFER_SIZE];
            snprintf(err, STR_BUFFER_SIZE, "dpotrf(M) failed (info: %d)", info);
            error_msg(err, 1);
        }

        // Wait for the previous gpu action to finish.
        sync_gpus(ngpus);

        START_SECTION("GPU_send_L", "Sending L to the GPUs");
        // GPU: Transfer L to the GPU(s) already.
        for(igpu = 0 ; igpu < ngpus ; ++igpu) {
            cudaSetDevice(igpu);
            // TODO: could do that in async.
            if((cu_status = cublasSetVector(L_gpu_bytes/sizeof(double), sizeof(double), M, 1, L_gpus[igpu], 1)) != CUBLAS_STATUS_SUCCESS) {
                char err[STR_BUFFER_SIZE];
                snprintf(err, STR_BUFFER_SIZE, "sending L to the GPU %d failed (info: %d)", igpu, cu_status);
                error_msg(err, 1);
            }
        }
        END_SECTION("GPU_send_L");

        /* XL := inv(L) XL */
        memcpy( XL, XL_orig, wXL * n * sizeof(double) );
        dtrsm_(LEFT, LOWER, NO_TRANS, NON_UNIT, &n, &wXL, &ONE, M, &n, XL, &n);

        START_SECTION2("WAIT_Y", "Waiting for current Y (%dth) of %d bytes (%d MB) from disk", j, aiocb_y_cur.aio_nbytes, aiocb_y_cur.aio_nbytes/1024.0/1024.0);
        /* Wait until current Y is available */
        /*printf("Waiting for %lu bytes read (Y)\n", aiocb_y_cur.aio_nbytes);*/
        fgls_aio_suspend( aiocb_y_cur_l, 1, NULL );
        END_SECTION("WAIT_Y");

        // Copy y for sigma2.score
        for( k = 0; k < n; k++ )
            scoreYmXB[k] = y_cur[k];

        /* y := inv(L) y */
        dtrsv_(LOWER, NO_TRANS, NON_UNIT, &n, M, &n, y_cur, &iONE);

        /* B_t := XL' * y */
        dgemv_(TRANS, &n, &wXL, &ONE, XL, &n, y_cur, &iONE, &ZERO, B_t, &iONE);

        /* V_tl := XL' * XL */
        dsyrk_(LOWER, TRANS, &wXL, &n, &ONE, XL, &n, &ZERO, V_tl, &wXL);
#ifdef VTRACE
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

        // At iteration iblock:
        // - The block iblock gets loaded from disk and sent to the GPU.
        // - The block iblock-1 gets computed on the GPU
        // - The block iblock-2 gets transferred back to and computed on the CPU
        int nblocks = m%x_b == 0 ? m/x_b : m/x_b+1;
        int iblock = 0;
        printf("going for %d blocks\n", nblocks);
        for (iblock = 0; iblock < nblocks+2; ++iblock)
        {
            printf("%d/%d\n", iblock, nblocks+1);
            // read i -> c
            // Read the next-but-one block of XR.
            // The last two iterations don't need to prefetch disk data anymore.
            if(iblock < nblocks) {
//                START_SECTION2("READ_X", "Starting the read of the Xr block (%d)", iblock);
                fgls_aio_read( &aiocb_x,
                               fileno( XR_fp ), X[hdd_cpu_buff],
                               MIN((size_t)x_b, m - ((size_t)x_b*iblock)) * wXR * n * sizeof(double),
                               ((off_t)x_b*iblock) * wXR * n * sizeof(double) );
//                END_SECTION("READ_X");
            }

            // complete alpha
            // Wait for the previous GPU computation (and transfer) to be done...
            // (The first GPU computation happens at i == 1 thus we first wait at i == 2)
            if(iblock >= 2) {
                START_SECTION2("GPU_trsm", "Waiting for the trsms (%d) on the GPUs to be done", iblock-2);
                for(igpu = 0 ; igpu < ngpus ; ++igpu) {
                    // cudaSetDevice(igpu); // Unnecessary according to "cuda_webinar_multi_gpu.pdf", p. 6
                    cudaStreamSynchronize(cu_comp_streams[igpu]);
                }
                END_SECTION("GPU_trsm");
            }

            // TODO: Would events be easier to read/understand? They still need multiple streams tho.

            if(iblock >= 1 && iblock <= nblocks) {
                // wait b -> beta
                // ...and wait for the sending of previous data to be done before...
                START_SECTION2("GPU_send_Xr", "Waiting for the sending of Xr (%d) parts to the GPUs to be done", iblock-1);
                for(igpu = 0 ; igpu < ngpus ; ++igpu) {
                    cudaStreamSynchronize(cu_trans_streams[igpu]);
                }
                END_SECTION("GPU_send_Xr");

                // compute beta
                // ...starting the next computation aswell as...
//                START_SECTION2("GPU_trsm", "Starting the trsms (%d) on the GPUs", iblock-1);
                /* XR := inv(L) XR */
                int curr_Xr_block_length = MIN(x_b, m - (iblock-1)*x_b);
                int rhss  = wXR * curr_Xr_block_length;

                for(igpu = 0 ; igpu < ngpus ; ++igpu) {
                    cudaSetDevice(igpu);
                    cublasSetStream(cu_handle, cu_comp_streams[igpu]);
                    if((cu_status = cublasDtrsm(cu_handle, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, n, rhss/ngpus, &ONE, L_gpus[igpu], n, Xr_beta_gpus[igpu], n)) != CUBLAS_STATUS_SUCCESS) {
                        char err[STR_BUFFER_SIZE];
                        snprintf(err, STR_BUFFER_SIZE, "computing the trsm on the GPU %d failed (info: %d)", igpu, cu_status);
                        error_msg(err, 1);
                    }
                }
//                END_SECTION("GPU_trsm");
            }

            size_t prev_Xr_block_length = MIN((size_t)x_b, m - ((size_t)x_b*(iblock-2)));

            // recv alpha -> a
            // ...the reception of the results of the previous computation.
            // (The first computation starts at i==1 thus the first result is available at i==2)
            if(iblock >= 2) {
//                START_SECTION2("GPU_recv_LXr", "Getting inv(L) * Xr (%d) back to the CPU", iblock-2);
                size_t prev_Xr_elems_per_device = prev_Xr_block_length*wXR*n/ngpus;
                for(igpu = 0 ; igpu < ngpus ; ++igpu) {
                    cudaSetDevice(igpu);
                    if((cu_status = cublasGetVectorAsync(prev_Xr_elems_per_device, sizeof(double), Xr_alpha_gpus[igpu], 1, X[cpu_gpu_buff] + igpu*prev_Xr_elems_per_device, 1, cu_trans_streams[igpu])) != CUBLAS_STATUS_SUCCESS) {
                        char err[STR_BUFFER_SIZE];
                        snprintf(err, STR_BUFFER_SIZE, "sending GPU %d's part of inv(L)*Xr to the CPU failed (info: %d)", igpu, cu_status);
                        error_msg(err, 1);
                    }
                }
//                END_SECTION("GPU_recv_LXr");
            }

            // wait c
            // Wait until current block of XR's is available
            // (In the last two iterations we don't need to read from disk anymore)
            if(iblock < nblocks) {
                START_SECTION2("WAIT_X", "Waiting for the next-but-one Xr block (%d) of %d bytes (%g MB)", iblock, aiocb_x.aio_nbytes, aiocb_x.aio_nbytes/1024.0/1024.0);
                fgls_aio_suspend( aiocb_x_l, 1, NULL );
                END_SECTION("WAIT_X");
            }

            // wait alpha -> a
            // ...and wait for the previous results to be sent back...
            if(iblock >= 2) {
                START_SECTION2("GPU_recv_LXr", "Waiting for getting inv(L)*Xr (%d) back to the CPU", iblock-2);
                for(igpu = 0 ; igpu < ngpus ; ++igpu) {
                    cudaStreamSynchronize(cu_trans_streams[igpu]);
                }
                END_SECTION("GPU_recv_LXr");

//                sync_write(X[cpu_gpu_buff], LXR_fp, prev_Xr_block_length*wXR*n, (size_t)x_b*(iblock-2)*wXR*n);
//                fflush(LXR_fp);
            }

            // send c -> alpha
            // ...so that now we can fill that GPU buffer with the next Xr block to compute
            if(iblock < nblocks) {
//                START_SECTION2("GPU_send_Xr", "Distributing the Xr block (%d) equally amongst the %d GPUs", iblock, ngpus);
                size_t Xr_elems_per_device = n*(wXR*x_b)/ngpus;
                for(igpu = 0 ; igpu < ngpus ; ++igpu) {
                    cudaSetDevice(igpu);
                    cublasSetStream(cu_handle, cu_trans_streams[igpu]);
                    if((cu_status = cublasSetVectorAsync(Xr_elems_per_device, sizeof(double), X[hdd_cpu_buff]+igpu*Xr_elems_per_device, 1, Xr_alpha_gpus[igpu], 1, cu_trans_streams[igpu])) != CUBLAS_STATUS_SUCCESS) {
                        char err[STR_BUFFER_SIZE];
                        snprintf(err, STR_BUFFER_SIZE, "sending part of Xr to the GPU %d failed (info: %d)", igpu, cu_status);
                        error_msg(err, 1);
                    }
                }
//                END_SECTION("GPU_send_Xr");
            }

            // CPU a -> a_hat

            if(iblock >= 2) {

                int x_inc = MIN(x_b, m - (iblock-2)*x_b);
                /* Set the number of threads for the multi-threaded BLAS to 1.
                 * The innermost loop is parallelized using OPENMP */
                #if defined GOTO
                setenv("GOTO_NUM_THREADS", "1", 1);
                #elif defined MKL
                setenv("MKL_NUM_THREADS",  "1", 1);
                #else
                setenv("OMP_NUM_THREADS",  "1", 1);
                #endif

                START_SECTION2("COMP_I", "Computing B and V (%d) on the CPU", iblock-2);
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
                            &ONE, &X[cpu_gpu_buff][i * wXR * n], &n, y_cur, &iONE, 
                            &ZERO, &Bij[wXL], &iONE);

                    /* Building V */
                    // Copy V_TL
                    for( k = 0; k < wXL; k++ )
                        dcopy_(&wXL, &V_tl[k*wXL], &iONE, &Vij[k*p], &iONE); // V_TL
                    // V_BL := XR' * XL
                    dgemm_("T", "N",
                            &wXR, &wXL, &n,
                            &ONE, &X[cpu_gpu_buff][i * wXR * n], &n, XL, &n,
                            &ZERO, &Vij[wXL], &p); // V_BL
                    // V_BR := XR' * XR
                    dsyrk_("L", "T", 
                            &wXR, &n, 
                            &ONE, &X[cpu_gpu_buff][i * wXR * n], &n, 
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
                END_SECTION("COMP_I");

                /* Wait until the previous blocks of B's and V's are written */
                if ( iblock > 2)
                {
                    START_SECTION2("WAIT_BV", "Waiting for the write of previous B and V (%d) to be done", iblock-3);
                    /*printf("Waiting for %lu bytes written (b)\n", aiocb_b_prev.aio_nbytes);*/
                    fgls_aio_suspend( aiocb_b_prev_l, 1, NULL );
                    /*printf("Waiting for %lu bytes written (v)\n", aiocb_v_prev.aio_nbytes);*/
                    fgls_aio_suspend( aiocb_v_prev_l, 1, NULL );
                    END_SECTION("WAIT_BV");
                }

                /* Write current blocks of B's and V's */
//                START_SECTION2("WRITE_BV", "Starting the writing of B and V (%d)", iblock-2);
                fgls_aio_write((struct aiocb*)aiocb_b_cur_l[0],
                                fileno( B_fp ), b_cur,
                                (size_t)x_inc * p * sizeof(double),
                                ((off_t)j * m * p + (off_t)x_b*(iblock-2) * p) * sizeof(double) );

                struct aiocb *aiocb_v_cur_p;
                aiocb_v_cur_p = (struct aiocb *) aiocb_v_cur_l[0];
                fgls_aio_write( aiocb_v_cur_p,
                                fileno( V_fp ), v_cur,
                                (size_t)x_inc * p * p * sizeof(double),
                                ((off_t)j * m * p * p + (off_t)x_b*(iblock-2) * p * p) * sizeof(double) );
//                END_SECTION("WRITE_BV");
                swap_aiocb( &aiocb_b_cur_l, &aiocb_b_prev_l );
                swap_buffers( &b_cur, &b_prev);
                swap_aiocb( &aiocb_v_cur_l, &aiocb_v_prev_l );
                swap_buffers( &v_cur, &v_prev);
            }

            /* Swap buffers */
            hdd_cpu_buff = (hdd_cpu_buff + 1) % 2;
            cpu_gpu_buff = (cpu_gpu_buff + 1) % 2;
            for(igpu = 0 ; igpu < ngpus ; ++igpu) {
                swap_buffers(Xr_alpha_gpus + igpu, Xr_beta_gpus + igpu);
            }
        }
        /* Swap buffers */
        swap_aiocb( &aiocb_y_cur_l, &aiocb_y_next_l );
        swap_buffers( &y_cur, &y_next);
    }

#ifdef VTRACE
    VT_USER_START("WAIT_ALL");
#endif
    /* Wait for the remaining IO operations issued */
      /*printf("Last Waiting");*/
    fgls_aio_suspend( aiocb_y_cur_l, 1, NULL );
    fgls_aio_suspend( aiocb_b_prev_l, 1, NULL );
    fgls_aio_suspend( aiocb_v_prev_l, 1, NULL );
#ifdef VTRACE
    VT_USER_END("WAIT_ALL");
#endif

    sync_gpus(ngpus); // Wait for the transfers of Xr to finish.

    START_SECTION("GPU_cleanup", "Waiting for all GPU stuff to end");

    // GPU: clean-up
    for(igpu = 0 ; igpu < ngpus ; ++igpu) {
        cudaSetDevice(igpu);
        cudaFree(L_gpus[igpu]);
        cudaFree(Xr_alpha_gpus[igpu]);
        cudaFree(Xr_beta_gpus[igpu]);
    }
    cublasDestroy(cu_handle);

    free(L_gpus);
    free(Xr_alpha_gpus);
    free(Xr_beta_gpus);

    end_section(&t_start, 0, "GPU_cleanup");

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

    cudaFreeHost(X[0]);
    cudaFreeHost(X[1]);
    cudaFreeHost(X[2]);
    free( Y[0] );
    free( Y[1] );
    free( B[0] );
    free( B[1] );
    free( V[0] );
    free( V[1] );

    free( aiocb_y_cur_l  );
    free( aiocb_y_next_l );
    free( aiocb_b_prev_l );
    free( aiocb_b_cur_l  );
    free( aiocb_v_prev_l );
    free( aiocb_v_cur_l  );

    return 0;
}

