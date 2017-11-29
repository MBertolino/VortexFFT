CC = gcc
LD = gcc
CFLAGS = -Wall -g -march=native #-Werror
LDFLAGS = -lfftw3 -lm 
INCLUDES =
RM = /bin/rm -f
OBJS = functions.o
PATCH = patch
FFT = swap

all:
	@echo "Usage: make patch"
patch: $(PATCH)
fft: $(FFT)


$(PATCH): patch.o $(OBJS)
	$(LD) -o $(PATCH) patch.o $(OBJS) $(LDFLAGS)

$(FFT): swap.o $(OBJS)
	$(LD) -o $(STATS) swap.o $(OBJS) $(LDFLAGS)

functions.o: functions.c functions.h
	$(CC) $(CFLAGS) $(INCLUDES) -c functions.c -std=c99

patch.o: patch.c
	$(CC) $(CFLAGS) $(INCLUDES) -c patch.c -std=c99

swap.o: swap.c
	$(CC) $(CFLAGS) $(INCLUDES) -c swap.c -std=c99


clean:
	$(RM) $(PATCH) $(FFT) *.o

