OBJDIR=bin
LDDIR=libs
SRCDIR=src
TESTDIR=test
INCDIR=include
CC=icc
CFLAGS=-I$(SRCDIR)/ -I$(TESTDIR)/ -I.

LIBS=$(LDDIR)/lapack_LINUX.a $(LDDIR)/libgoto2_penrynp-r1.07.a
#LIBS=$(LDDIR)/lapack_LINUX.a $(LDDIR)/libgoto2.a
#LIBS=$(FLAGS_MKL_LINKER)

# all src files which do not contain a 'main' definition'
NON_MAIN_SRC := src/eigenDec.c \
		src/fgls.c \
		src/io.c \
		src/fgls_eigen.c \
		src/mod_x_y.c \
		test/test_framework.c

NON_MAIN_OBJ := \
	$(patsubst %.c,$(OBJDIR)/%.o,$(NON_MAIN_SRC))

.all: main read_test read_blocksize_test read_tests write_test write_blocksize_test write_tests write_h_file clean
.PHONY: .all 

all: main


$(OBJDIR)/%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)

$(OBJDIR)/driver.x: bin/src/driver.o $(NON_MAIN_OBJ)
	$(CC) $^ $(LIBS)  -lgfortran -lpthread  $(CFLAGS) -o $@

$(OBJDIR)/read_test.x: bin/test/read_test.o $(NON_MAIN_OBJ)
	$(CC) $^ $(LIBS)  -lgfortran -lpthread  $(CFLAGS) -o $@

$(OBJDIR)/read_blocksize_test.x: bin/test/read_blocksize_test.o $(NON_MAIN_OBJ)
	$(CC) $^ $(LIBS)  -lgfortran -lpthread  $(CFLAGS) -o $@

$(OBJDIR)/write_test.x: bin/test/write_test.o $(NON_MAIN_OBJ)
	$(CC) $^ $(LIBS)  -lgfortran -lpthread  $(CFLAGS) -o $@

$(OBJDIR)/write_blocksize_test.x: bin/test/write_blocksize_test.o $(NON_MAIN_OBJ)
	$(CC) $^ $(LIBS)  -lgfortran -lpthread  $(CFLAGS) -o $@

$(OBJDIR)/write_h_file.x: bin/src/write_h_file.o
	$(CC) $^ $(CFLAGS) -o $@


main: $(OBJDIR)/driver.x 
read_test: $(OBJDIR)/read_test.x 
read_blocksize_test: $(OBJDIR)/read_blocksize_test.x 
read_tests: read_test read_blocksize_test 
write_test: $(OBJDIR)/write_test.x 
write_blocksize_test: $(OBJDIR)/write_blocksize_test.x 
write_tests: write_test write_blocksize_test 
write_h_file: $(OBJDIR)/write_h_file.x

clean:
	rm -rf $(OBJDIR); 
	rm -f $(SRCDIR)/*~ $(TESTDIR)/*~ core;


