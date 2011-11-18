CC = vtcc
#CFLAGS = -D_FILE_OFFSET_BITS=64 -DVTRACE -vt:cc gcc -vt:mt -vt:inst manual -pthread -fopenmp -DMKL
#CFLAGS =  -pthread -fopenmp -DMKL -DVTRACE -vt:cc gcc -vt:mt -vt:inst manual -vt:noopari -vt:verbose
#CFLAGS = -D_FILE_OFFSET_BITS=64 -pthread -openmp -DMKL -DVTRACE -vt:cc icc -vt:mt -vt:inst manual -vt:verbose
CFLAGS = -D_FILE_OFFSET_BITS=64 -pthread -DMKL -DVTRACE -vt:cc gcc -vt:mt -vt:inst manual -vt:verbose

LDFLAGS = -L$(MKLROOT)/lib/intel64/ -L$(INTELROOT)/compiler/lib/intel64
LDLIBS = -lrwthmkl -liomp5 -lifcoremt
