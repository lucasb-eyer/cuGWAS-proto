include ./make.inc

SRCDIR = ./src
DRIVER = ./bin/driver.x

CFLAGS+=-g -Wall -I$(SRCDIR)/
LDLIBS += -lm -lrt

SRCS = $(SRCDIR)/common.c $(SRCDIR)/driver.c $(SRCDIR)/fgls_chol.c $(SRCDIR)/fgls_eigen.c $(SRCDIR)/io.c $(SRCDIR)/timing.c
OBJS = $(SRCS:.c=.o)
EXE=./bin/driver.x

.PHONY: all clean 

all: $(EXE)

$(EXE): $(OBJS)
	$(CC) $(CFLAGS) $^ $(LDFLAGS) $(LDLIBS) -o $@

gpu: $(OBJS) $(SRCDIR)/fgls_chol_gpu.o
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -L$(CUDA_ROOT)/lib64 $(LDLIBS) -lcublas -o $(EXE)

$(SRCDIR)/fgls_chol_gpu.o: $(SRCDIR)/fgls_chol_gpu.c
	$(CC) $(CFLAGS) -I$(CUDA_ROOT)/include -DFGLS_WITH_GPU -DFGLS_GPU_SERIAL -c $^ -o $@
$(SRCDIR)/driver.o: $(SRCDIR)/driver.c
	$(CC) $(CFLAGS) -I$(CUDA_ROOT)/include -DFGLS_WITH_GPU -c $^ -o $@

clean:
	$(RM) $(OBJS) $(SRCDIR)/fgls_chol_gpu.o;
	$(RM) $(EXE);
	$(RM) $(SRCDIR)/*mod*
	$(RM) $(SRCDIR)/*opari*
