#include "fgls.h"
#include "io.h"

#include <stdio.h>
#include <stdlib.h>
/*#include <unistd.h>*/
#include <math.h>

#if TIMING
#include <sys/time.h>
#include <time.h>
#endif // TIMING

double compare(double* a, double* b, int size) {
  double out = 0.0;
  int i;
  for (i = 0; i < size; i++) {
    if((i % 1000) == 0 && fabs(a[i] - b[i]) > 1e-13)
		printf("Difference at %d: %e\n", i, fabs(a[i] - b[i]));
    if(out < fabs(a[i] - b[i]))
      out = fabs(a[i] - b[i]);
  }
  return out;
}

#if TIMING
long get_diff_ms(struct timeval *start_time, struct timeval *end_time) {
  long seconds = end_time->tv_sec - start_time->tv_sec;
  long useconds = end_time->tv_usec - start_time->tv_usec;
  return ((seconds) * 1000 + useconds/1000.0) +0.5;
}
long get_diff_us(struct timeval *start_time, struct timeval *end_time) {
  long seconds = end_time->tv_sec - start_time->tv_sec;
  long useconds = end_time->tv_usec - start_time->tv_usec;
  return ((seconds) * 1000000 + useconds);
}
long get_diff_ns(struct timespec *start, struct timespec *end)
{
	/*timespec temp;*/
	long seconds, nseconds;
	long nsecs_in_sec = 1e9;
	if ( (end->tv_nsec - start->tv_nsec) < 0 ) {
		/*temp.tv_sec = end.tv_sec-start.tv_sec-1;*/
		/*temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;*/
        seconds = end->tv_sec - start->tv_sec - 1;
		nseconds = nsecs_in_sec - (end->tv_nsec - start->tv_nsec);
		return nsecs_in_sec * seconds + nseconds;
	} else {
		/*temp.tv_sec = end.tv_sec-start.tv_sec;*/
		/*temp.tv_nsec = end.tv_nsec-start.tv_nsec;*/
        seconds = end->tv_sec - start->tv_sec;
		nseconds = end->tv_nsec - start->tv_nsec;
		return nsecs_in_sec * seconds + nseconds;
	}
}
#endif // TIMING

void swap_buffers(double** b1, double** b2) {
  double* temp;
  temp = *b1;
  *b1 = *b2;
  *b2 = temp;
}

/*void read_x(double* buf, int index, const problem_args* args) {*/
/*if(!x_file) {*/
/*printf("x_file not initialized. Exiting...\n");*/
/*exit(-1);*/
/*}*/
/*int x_inc = MIN(args->x_b, args->m - args->x_b*index);*/
/*read(buf, x_file, args->p*args->n*x_inc, index*args->x_b*args->p*args->n);*/
/*}*/

/*void read_phi(double* buf, int index, const problem_args* args) {*/
/*if(!phi_file) {*/
/*printf("phi_file not initialized. Exiting...\n");*/
/*exit(-1);*/
/*}*/
/*read(buf, phi_file, args->n*args->n, 0);*/
/*}*/

/*void read_y(double* buf, int index, const problem_args* args) {*/
/*if(!y_file) {*/
/*printf("y_file not initialized. Exiting...\n");*/
/*exit(-1);*/
/*}*/
/*int y_inc = MIN(args->y_b, args->t - args->y_b*index); */
/*read(buf, y_file, args->n*y_inc, index*args->y_b*args->n);*/
/*}*/

/*void read_h(double* buf, int index, const problem_args* args) {*/
/*if(!h_file) {*/
/*printf("h_file not initialized. Exiting...\n");*/
/*exit(-1);*/
/*}*/
/*read(buf, h_file, args->t, 0);*/
/*}*/

/*int return_buffer_index(double** buffers, int size, double* cur) {*/
/*int i;*/
/*for(i = 0; i < size; i++) {*/
/*if (buffers[i] == cur) {*/
/*return i;*/
/*}*/
/*}*/
/*return -1;*/
/*}*/


/*void write_b(double* buf, int r, int s, const problem_args* args) {*/
/*if(!b_file) {*/
/*printf("b_file not initialized. Exiting...\n");*/
/*exit(-1);*/
/*}*/
/*int y_inc, x_inc, j, buffer_index, file_index;*/
/*y_inc = MIN(args->y_b, args->t - args->y_b*s);*/
/*x_inc = MIN(args->x_b, args->m - args->x_b*r);*/

/*for (j = 0; j < y_inc; j++) {*/
/*buffer_index = x_inc*args->p*j;*/
/*file_index = args->m*args->p*(args->y_b*s + j) + r*args->x_b*args->p;*/
/*write(&buf[buffer_index], b_file, args->p*x_inc, file_index);*/
/*}*/
/*}*/

/*void write_x(double* buf, int r, const problem_args* args) {*/
/*if(!x_tmp_file) {*/
/*printf("x_tmp_file not initialized. Exiting...\n");*/
/*exit(-1);*/
/*}*/
/*int x_inc, i, buffer_index, file_index;*/
/*x_inc = MIN(args->x_b, args->m - args->x_b*r);*/

/*file_index = args->p*args->n*r*args->x_b;*/
/*write(buf, x_tmp_file, args->p*args->n*x_inc, file_index);*/
/*}*/

/*void write_y(double* buf, int s, const problem_args* args) {*/
/*if(!y_tmp_file) {*/
/*printf("y_tmp_file not initialized. Exiting...\n");*/
/*exit(-1);*/
/*}*/
/*int y_inc, j, buffer_index, file_index;*/
/*y_inc = MIN(args->y_b, args->t - args->y_b*s);*/

/*file_index = args->n*s*args->y_b;*/
/*write(buf, y_tmp_file, args->n*y_inc, file_index);*/
/*}*/

