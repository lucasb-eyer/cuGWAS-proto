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
C_FILES_SRC := $(wildcard $(SRCDIR)/*.c)
C_FILES_TEST := $(wildcard $(TESTDIR)/*.c)
OBJ_FILES_SRC := \
	$(patsubst $(SRCDIR)/%.c,$(OBJDIR)/$(SRCDIR)/%.o,$(C_FILES_SRC))
OBJ_FILES_TEST := \
	$(patsubst $(TESTDIR)/%.c,$(OBJDIR)/$(TESTDIR)/%.o,$(C_FILES_TEST))

all: $(OBJDIR)/driver.x $(OBJ_FILES_TEST) $(OBJ_FILES_SRC)

$(OBJDIR)/$(SRCDIR)/%.o: $(SRCDIR)/%.c $(DEPS)
	mkdir -p $(OBJDIR)/$(SRCDIR);
	$(CC) -c -o $@ $< $(CFLAGS)
$(OBJDIR)/$(TESTDIR)/%.o: $(TESTDIR)/%.c $(DEPS)
	mkdir -p $(OBJDIR)/$(TESTDIR);
	$(CC) -c -o $@ $< $(CFLAGS)

$(OBJDIR)/driver.x: $(OBJ_FILES_TEST) $(OBJ_FILES_SRC)
	$(CC) $^ $(LIBS)  -lgfortran -lpthread  $(CFLAGS) -o $@
#	$(CC) $^ $(LIBS)  $(CFLAGS) -o $@

clean:
	rm -rf $(OBJDIR); 
	rm -f $(SRCDIR)/*~ $(TESTDIR)/*~ core;
