ODIR=bin
LDIR=lib
SDIR=src
IDIR=$(SDIR)

CC=gcc
CFLAGS=-Wall -g -I$(IDIR) 

LIBS=-lpthread
C_FILES := $(wildcard $(SDIR)/*.c)
OBJ_FILES := $(patsubst $(SDIR)/%.c,$(ODIR)/%.o,$(C_FILES))

_DEPS = fgls.h io.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

all: $(ODIR)/driver.x $(OBJ_FILES)

$(ODIR)/%.o: $(SDIR)/%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(ODIR)/driver.x: $(OBJ_FILES)
	gcc -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm -f $(ODIR)/*.o $(SDIR)/*~ core
