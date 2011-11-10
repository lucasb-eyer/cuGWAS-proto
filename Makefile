include ./make.inc

SRCDIR = ./src
DRIVER = ./bin/driver.x

CFLAGS+=-g -Wall -I$(SRCDIR)/
LDLIBS += -lm -lrt

SRCS = $(wildcard $(SRCDIR)/*.c)
OBJS = $(SRCS:.c=.o)
EXE=./bin/driver.x


.PHONY: all clean 

all: $(EXE)

$(EXE): $(OBJS)
	$(CC) $(CFLAGS) $^ $(LDFLAGS) $(LDLIBS) -o $@

clean:
	$(RM) $(OBJS); 
	$(RM) $(EXE); 
	$(RM) $(SRCDIR)/*mod*
	$(RM) $(SRCDIR)/*opari*
