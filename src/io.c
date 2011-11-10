#include <stdio.h>
#include <stdlib.h>
#include <strings.h>

#include "io.h"

void sync_read(double *buffer, FILE* fp, size_t size, size_t start) 
{
	size_t out;

	fseek(fp, start*sizeof(double), SEEK_SET);
	out = fread(buffer, sizeof(double), size, fp);
	if ( out != size )
	{
		printf("[ERROR] "__FILE__ ": Read %d elements - %d expected\n", (int)out, (int)size);
	}

	return; 
}

void sync_write(double* buffer, FILE* fp, size_t size, size_t start) 
{
	size_t out;

	fseek(fp, start * sizeof(double), SEEK_SET);
	out = fwrite(buffer, sizeof(double), size, fp);
	if ( out != size ) 
	{
		printf("[ERROR] "__FILE__ ": Wrote %d elements - %d expected\n", (int)out, (int)size);
	}

	return; 
}

void fgls_aio_read(struct aiocb *aiocb, int fildes, void *buf, size_t nbytes, off_t offset)
{
	bzero( (char *)aiocb,  sizeof(struct aiocb) );

	/*printf("Reading %zu bytes\n", nbytes);*/
	/*printf("Reading offset: %jd\n", offset);*/
	/*fflush(stdout);*/
	aiocb->aio_fildes = fildes;
	aiocb->aio_buf = buf;
	aiocb->aio_nbytes = nbytes;
	aiocb->aio_offset = offset;
	/*aiocb_y_cur.aio_flag = AIO_RAW;*/

	if ( aio_read( aiocb ) != 0 )
	{
		perror("aio_read error");
		exit( EXIT_FAILURE );
	}

	return;
}

void fgls_aio_write(struct aiocb *aiocb, int fildes, void *buf, size_t nbytes, off_t offset)
{
	bzero( (char *)aiocb,  sizeof(struct aiocb) );

	aiocb->aio_fildes = fildes;
	aiocb->aio_buf = buf;
	aiocb->aio_nbytes = nbytes;
	aiocb->aio_offset = offset;
	/*aiocb_y_cur.aio_flag = AIO_RAW;*/

	if ( aio_write( aiocb ) != 0 )
	{
		perror("aio_write error");
		exit( EXIT_FAILURE );
	}

	return;
}

void fgls_aio_suspend(const struct aiocb * const cblist[], 
		              int n, const struct timespec *timeout)
{
	if ( aio_suspend( cblist, n, timeout ) != 0 )
	{
		perror("aio_suspend error");
		exit( EXIT_FAILURE );
	}
	/*printf( "%zu\n", cblist[0]->aio_nbytes );*/
	/*fflush(stdout);*/
	if ( aio_return( (struct aiocb *)cblist[0] ) != cblist[0]->aio_nbytes )
	{
		perror("aio_suspend error - data not read completely");
		exit( EXIT_FAILURE );
	}
}
