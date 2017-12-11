CC = gcc
LD = gcc
CFLAGS = -Wall -O3 -march=native #-Werror
LDFLAGS = -lfftw3 -lm 
INCLUDES =
RM = /bin/rm -f
OBJS = orig_functions.o fft_functions.o misc.o
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

orig_functions.o: orig_functions.c orig_functions.h
	$(CC) $(CFLAGS) $(INCLUDES) -c orig_functions.c -std=c99

fft_functions.o: fft_functions.c fft_functions.h
	$(CC) $(CFLAGS) $(INCLUDES) -c fft_functions.c -std=c99

misc.o: misc.c misc.h
	$(CC) $(CFLAGS) $(INCLUDES) -c misc.c -std=c99

patch.o: patch.c
	$(CC) $(CFLAGS) $(INCLUDES) -c patch.c -std=c99


cleanres:
	$(RM) ../results/*.txt

clean:
	$(RM) $(PATCH) $(FFT) *.o

