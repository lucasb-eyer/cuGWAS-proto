OBJDIR=bin
LDDIR=lib
SRCDIR=src
TESTDIR=test
INCDIR=$(SDIR)

CC=gcc
CFLAGS=-Wall -g -I$(SRCDIR) -I$(TESTDIR) 

LIBS=-lpthread
C_FILES := $(wildcard $(SRCDIR)/*.c)
C_FILES := $(wildcard $(TESTDIR)/*.c)
OBJ_FILES := $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(C_FILES))
OBJ_FILES := $(patsubst $(TESTDIR)/%.c,$(OBJDIR)/%.o,$(C_FILES))

_DEPS = fgls.h io.h test_framework.h read_test.h
DEPS := $(patsubst %,$(SRCDIR)/%,$(_DEPS))
DEPS := $(patsubst %,$(TESTDIR)/%,$(_DEPS))

all: $(OBJDIR)/driver.x $(OBJ_FILES)

$(OBJDIR)/%.o: $(SRCDIR)/%.c $(TESTDIR)/%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(OBJDIR)/driver.x: $(OBJ_FILES)
	gcc -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm -f $(OBJDIR)/*.o $(SRCDIR)/*~ $(TESTDIR)/*~ core
