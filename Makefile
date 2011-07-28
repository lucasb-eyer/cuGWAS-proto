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
NON_MAIN_SRC := src/bio_eigen.c \
		src/bio_chol.c \
		src/eigenDec.c \
		src/fgls.c \
		src/io.c \
		src/m_traversal_chol.c \
		src/m_traversal_eigen.c \
		src/mod_x_y.c \
		src/t_traversal_chol.c \
		src/t_traversal_eigen.c \
		test/test_framework.c

NON_MAIN_OBJ := \
	$(patsubst %.c,$(OBJDIR)/%.o,$(NON_MAIN_SRC))

all: main
.all: main read_test clean
.PHONY: .all

$(OBJDIR)/%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)

$(OBJDIR)/driver.x: bin/src/driver.o $(NON_MAIN_OBJ)
	$(CC) $^ $(LIBS)  -lgfortran -lpthread  $(CFLAGS) -o $@
#	$(CC) $^ $(LIBS)  $(CFLAGS) -o $@

$(OBJDIR)/read_test.x: bin/test/read_test.o $(NON_MAIN_OBJ)
	$(CC) $^ $(LIBS)  -lgfortran -lpthread  $(CFLAGS) -o $@

main: $(OBJDIR)/driver.x
read_test: $(OBJDIR)/read_test.x
clean:
	rm -rf $(OBJDIR); 
	rm -f $(SRCDIR)/*~ $(TESTDIR)/*~ core;


