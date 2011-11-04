#ifndef IO_H
#define IO_H

#include <stdio.h>
#include <aio.h>

void sync_write(double *buffer, FILE* file, size_t size, size_t start);

void sync_read(double *buffer, FILE* file, size_t size, size_t start);

void fgls_aio_read(struct aiocb *aiocb, int fildes, void *buf, size_t nbytes, off_t offset);

void fgls_aio_write(struct aiocb *aiocb, int fildes, void *buf, size_t nbytes, off_t offset);

void fgls_aio_suspend(const struct aiocb * const cblist[], 
		              int n, const struct timespec *timeout);

#endif
