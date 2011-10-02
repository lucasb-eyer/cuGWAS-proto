OBJDIR=bin
LDDIR=libs
SRCDIR=src
TESTDIR=test
INCDIR=include
CC=gcc
#CC=icc
CFLAGS= -g -Wall -pthread -I$(SRCDIR)/ -I$(TESTDIR)/ -I.

#CC=vtcc
#CFLAGS+= -g -DVTRACE -vt:mt -vt:inst manual -pthread -vt:verbose

#LIBS=$(LDDIR)/lapack_LINUX.a $(LDDIR)/libgoto2_penrynp-r1.07.a
#LIBS=$(LDDIR)/lapack_LINUX.a $(LDDIR)/libgoto2.a
#LIBS=$(FLAGS_MKL_LINKER)
LIBS = -lm $(HOME)/libs/lapack-3.3.1/lapack_LINUX.a $(HOME)/libs/GotoBLAS2/libgoto2.a -lgfortran #-lpthread
#LIBS = -ldl -lm $(HOME)/libs/lapack-3.3.1/lapack_LINUX.a $(HOME)/libs/GotoBLAS2/libgoto2.a -lgfortran -lpthread
#LIBPATHS=$(HOME)/libs/GotoBLAS2/

# all src files which do not contain a 'main' definition'
NON_MAIN_SRC := src/timing.c \
		src/io.c \
		src/fgls_eigen.c \
		src/preloop.c \
		test/test_framework.c

NON_MAIN_OBJ := \
	$(patsubst %.c,$(OBJDIR)/%.o,$(NON_MAIN_SRC))

.all: main read_test read_blocksize_test read_tests write_test write_blocksize_test write_tests write_h_file clean
.PHONY: .all 

all: main


$(OBJDIR)/%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)

$(OBJDIR)/driver.x: bin/src/driver.o $(NON_MAIN_OBJ)
	$(CC) $^ $(LIBS) -L$(LIBPATHS) $(CFLAGS) -lrt -o $@

$(OBJDIR)/read_test.x: bin/test/read_test.o $(NON_MAIN_OBJ)
	$(CC) $^ $(LIBS) -L$(LIBPATHS) $(CFLAGS) -o $@

$(OBJDIR)/read_blocksize_test.x: bin/test/read_blocksize_test.o $(NON_MAIN_OBJ)
	$(CC) $^ $(LIBS) -L$(LIBPATHS) $(CFLAGS) -o $@

$(OBJDIR)/write_test.x: bin/test/write_test.o $(NON_MAIN_OBJ)
	$(CC) $^ $(LIBS) -L$(LIBPATHS) $(CFLAGS) -o $@

$(OBJDIR)/write_blocksize_test.x: bin/test/write_blocksize_test.o $(NON_MAIN_OBJ)
	$(CC) $^ $(LIBS) -L$(LIBPATHS)  $(CFLAGS) -o $@

$(OBJDIR)/write_h_file.x: bin/src/write_h_file.o
	$(CC) $^ $(CFLAGS) -o $@


main: dirs $(OBJDIR)/driver.x
dirs:
	mkdir -p bin/src/
	mkdir -p bin/test/
read_test: $(OBJDIR)/read_test.x 
read_blocksize_test: $(OBJDIR)/read_blocksize_test.x 
read_tests: read_test read_blocksize_test 
write_test: $(OBJDIR)/write_test.x 
write_blocksize_test: $(OBJDIR)/write_blocksize_test.x 
write_tests: write_test write_blocksize_test 
write_h_file: $(OBJDIR)/write_h_file.x

bin/src/driver.o: src/driver.c src/io.h src/timing.h src/fgls_eigen.h
bin/src/fgls_eigen.o: src/fgls_eigen.c src/blas.h src/lapack.h src/io.h \
 src/timing.h src/fgls_eigen.h
bin/src/io.o: src/io.c src/io.h
bin/src/preloop.o: src/preloop.c src/options.h src/blas.h src/lapack.h \
 src/timing.h src/io.h
bin/src/timing.o: src/timing.c src/io.h src/timing.h

clean:
	rm -rf $(OBJDIR); 
	rm -f $(SRCDIR)/*~ $(TESTDIR)/*~ core;


