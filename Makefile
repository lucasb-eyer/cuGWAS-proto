
# indicate where the object files are to be created
CC         := gcc 
LINKER     := $(CC) -lpthread
#CFLAGS     := -O3 -Wall
CFLAGS     := -Wall -g

HEADERS := io.h

TEST_OBJS  := main.o \
	      io.o	

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@


driver.x: $(TEST_OBJS) $(HEADERS)
	$(LINKER) $(TEST_OBJS) -o main.x
test: driver.x
	./main.x a b c < 1 1 1 1 4 4
clean:
	rm -f *.o *~ core *.x

