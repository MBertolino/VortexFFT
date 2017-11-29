CC = gcc
LD = gcc
CFLAGS = -Wall -g -march=native #-Werror
LDFLAGS = -lm 
INCLUDES =
RM = /bin/rm -f
OBJS = contiguous_functions.o
PATCH = patch
FFT = swap

all:
	@echo "Usage: make patch"
patch: $(PATCH)
fft: $(FFT)


$(PATCH): contiguous_patch.o $(OBJS)
	$(LD) -o $(PATCH) contiguous_patch.o $(OBJS) $(LDFLAGS)

$(FFT): swap.o $(OBJS)
	$(LD) -o $(STATS) swap.o $(OBJS) $(LDFLAGS)

contiguous_functions.o: contiguous_functions.c contiguous_functions.h
	$(CC) $(CFLAGS) $(INCLUDES) -c contiguous_functions.c -std=c99

contiguous_patch.o: contiguous_patch.c
	$(CC) $(CFLAGS) $(INCLUDES) -c contiguous_patch.c -std=c99

swap.o: swap.c
	$(CC) $(CFLAGS) $(INCLUDES) -c swap.c -std=c99


clean:
	$(RM) $(PATCH) $(FFT) *.o

