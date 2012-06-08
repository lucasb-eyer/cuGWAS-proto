#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/time.h>

#include "common.h"

/*void error_msg(char *file, int line, char *msg, int abort)*/
void error_msg(char *msg, int abort)
{
    /*fprintf(stderr, "[Error] %s(line %d): %s\n", file, line, msg);*/
    fprintf(stderr, "[Error]: %s\n", msg);
    if (abort)
        exit(EXIT_FAILURE);
}

void * fgls_malloc_impl( const char* file, long line, size_t size )
{
    void *buf;

    if ( (buf = malloc( size )) == NULL ) {
        fprintf(stderr, "Couldn't allocate %ld bytes of memory in %s:%ld\n", size, file, line);
        exit(EXIT_FAILURE);
    }

    memset(buf, 0, size);

    return buf;
}

int read_clock(struct timeval *t)
{
    return gettimeofday(t, NULL);
}

int elapsed_time(struct timeval *start, struct timeval *end)
{
    return (int)(end->tv_sec - start->tv_sec)*1e6 + (int)(end->tv_usec - start->tv_usec);
}

void sleep_seconds(double seconds)
{
    double real = 0;
    struct timespec ts;
    ts.tv_nsec = (long)(modf(seconds, &real) * 1000000.0);
    ts.tv_sec = (time_t)real;
    nanosleep(&ts, 0);
}

void initialize_config(
        FGLS_config_t *cf,
        int n, int p, int m, int t, int wXL, int wXR,
        int x_b, int y_b, int num_threads,
        char *Phi_path, char *h2_path, char *sigma2_path,
        char *XL_path, char *XR_path, char *Y_path,
        char *B_path, char *V_path
)
{
    // Problem dimensions
    cf->n = n;
    cf->p = p;
    cf->m = m;
    cf->t = t;
    cf->wXR = wXR;
    cf->wXL = wXL;
    // Algorithm parameters
    cf->x_b = x_b;
    cf->y_b = y_b;
    cf->NUM_COMPUTE_THREADS = num_threads;
    // In/Out Files
    /*strncpy( cf->Phi_path,    Phi_path,    STR_BUFFER_SIZE );*/
    /*strncpy( cf->h_path,      h2_path,     STR_BUFFER_SIZE );*/
    /*strncpy( cf->sigma_path , sigma2_path, STR_BUFFER_SIZE );*/
    /*strncpy( cf->XL_path,     XL_path,     STR_BUFFER_SIZE );*/
    /*strncpy( cf->XR_path,     XR_path,     STR_BUFFER_SIZE );*/
    /*strncpy( cf->Y_path,      Y_path,      STR_BUFFER_SIZE );*/
    /*strncpy( cf->B_path,      B_path,      STR_BUFFER_SIZE );*/
    /*strncpy( cf->V_path,      V_path,      STR_BUFFER_SIZE );*/
    sprintf( cf->Phi_path,    "%s", Phi_path );
    sprintf( cf->h_path,      "%s", h2_path );
    sprintf( cf->sigma_path , "%s", sigma2_path );
    sprintf( cf->XL_path,     "%s", XL_path );
    sprintf( cf->XR_path,     "%s", XR_path );
    sprintf( cf->Y_path,      "%s", Y_path );
    sprintf( cf->B_path,      "%s", B_path );
    sprintf( cf->V_path,      "%s", V_path );
    /*strcpy( cf->Phi_path,    Phi_path );*/
    /*strcpy( cf->h_path,      h2_path );*/
    /*strcpy( cf->sigma_path , sigma2_path );*/
    /*strcpy( cf->XL_path,     XL_path );*/
    /*strcpy( cf->XR_path,     XR_path );*/
    /*strcpy( cf->Y_path,      Y_path );*/
    /*strcpy( cf->B_path,      B_path );*/
    /*strcpy( cf->V_path,      V_path );*/
    // Temporary files
    snprintf( cf->ZtXL_path, STR_BUFFER_SIZE, "%s.tmp", XL_path );
    snprintf( cf->ZtXR_path, STR_BUFFER_SIZE, "%s.tmp", XR_path );
    snprintf( cf->ZtY_path,  STR_BUFFER_SIZE, "%s.tmp", Y_path );

    return;
}

void swap_buffers(double** b1, double** b2) 
{
    double* tmp;

    tmp = *b1;
    *b1 = *b2;
    *b2 = tmp;
}

void swap_aiocb(const struct aiocb ***x, const struct aiocb ***y)
{
    const struct aiocb **tmp = *x;
    *x = *y;
    *y = tmp;
}

