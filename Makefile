OBJDIR=bin
LDDIR=lib
SRCDIR=src
TESTDIR=test
INCDIR=$(SDIR)

CC=gcc
CFLAGS=-Wall -g -I. -I$(SRCDIR)/ -I$(TESTDIR)/ 

LIBS=-lpthread
C_FILES_SRC := $(wildcard $(SRCDIR)/*.c)
C_FILES_TEST := $(wildcard $(TESTDIR)/*.c)
OBJ_FILES_SRC := \
	$(patsubst $(SRCDIR)/%.c,$(OBJDIR)/$(SRCDIR)/%.o,$(C_FILES_SRC))
OBJ_FILES_TEST := \
	$(patsubst $(TESTDIR)/%.c,$(OBJDIR)/$(TESTDIR)/%.o,$(C_FILES_TEST))

DEPS = $(SRCDIR)/fgls.h $(SRCDIR)/io.h \
	$(TESTDIR)/test_framework.h 

all: $(OBJDIR)/driver.x $(OBJ_FILES_TEST) $(OBJ_FILES_SRC)

$(OBJDIR)/$(SRCDIR)/%.o: $(SRCDIR)/%.c $(DEPS)
	mkdir -p $(OBJDIR)/$(SRCDIR);
	$(CC) -c -o $@ $< $(CFLAGS)
$(OBJDIR)/$(TESTDIR)/%.o: $(TESTDIR)/%.c $(DEPS)
	mkdir -p $(OBJDIR)/$(TESTDIR);
	$(CC) -c -o $@ $< $(CFLAGS)

$(OBJDIR)/driver.x: $(OBJ_FILES_TEST) $(OBJ_FILES_SRC)
	gcc -o $@ $^ $(CFLAGS) $(LIBS)

echo:
	echo $(C_FILES);
	echo $(OBJ_FILES);
	echo $(DEPS);
clean:
	rm -rf $(OBJDIR); 
	rm -f $(SRCDIR)/*~ $(TESTDIR)/*~ core;
