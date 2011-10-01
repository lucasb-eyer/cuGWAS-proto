#include <stdio.h>

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

