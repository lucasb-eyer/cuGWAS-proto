CC = vtcc
#CFLAGS = -D_FILE_OFFSET_BITS=64 -DVTRACE -vt:cc gcc -vt:mt -vt:inst manual -pthread -fopenmp -DMKL
#CFLAGS =  -pthread -fopenmp -DMKL -DVTRACE -vt:cc gcc -vt:mt -vt:inst manual -vt:noopari -vt:verbose
CFLAGS =  -pthread -fopenmp -DMKL -DVTRACE -vt:cc gcc -vt:mt -vt:inst manual -vt:verbose

LDFLAGS = -L$(MKLROOT)/lib/ -L/opt/intel/composerxe-2011.1.122/compiler/lib/
LDLIBS = -lmkl_intel -lmkl_intel_thread -lmkl_core -liomp5 -lmkl_lapack95_ilp64 # -lguide 
