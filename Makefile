include ./make.inc

SRCDIR = ./src
DRIVER = ./bin/driver.x

CFLAGS+=-g -Wall -I$(SRCDIR)/
LDLIBS += -lm -lrt

SRCS = $(SRCDIR)/common.c $(SRCDIR)/driver.c $(SRCDIR)/fgls_chol.c $(SRCDIR)/fgls_eigen.c $(SRCDIR)/io.c $(SRCDIR)/timing.c
OBJS = $(SRCS:.c=.o)
EXE=./bin/driver.x
EXE_GPU=./bin/driver.gpu.x

.PHONY: all clean

all: $(EXE)

$(EXE): $(OBJS)
	$(CC) $(CFLAGS) -I$(CUDA_ROOT)/include -c $(SRCDIR)/driver.c -o $(SRCDIR)/driver.o
	$(CC) $(CFLAGS) $^ $(LDFLAGS) $(LDLIBS) -o $@

gpu: $(OBJS) $(SRCDIR)/fgls_chol_gpu.o
	$(CC) $(CFLAGS) -I$(CUDA_ROOT)/include -DFGLS_WITH_GPU -DFGLS_GPU_SERIAL -c $(SRCDIR)/fgls_chol_gpu.c -o $(SRCDIR)/fgls_chol_gpu.o
	$(CC) $(CFLAGS) -I$(CUDA_ROOT)/include -DFGLS_WITH_GPU -c $(SRCDIR)/driver.c -o $(SRCDIR)/driver.o
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -L$(CUDA_ROOT)/lib64 $(LDLIBS) -lcublas -o $(EXE_GPU)

clean:
	$(RM) $(OBJS) $(SRCDIR)/fgls_chol_gpu.o;
	$(RM) $(EXE) $(EXE_GPU);
	$(RM) $(SRCDIR)/*mod*
	$(RM) $(SRCDIR)/*opari_GPU*
