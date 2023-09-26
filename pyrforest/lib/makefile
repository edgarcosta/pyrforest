CC = gcc
##### using -fPIC slows things down by a few percent, not a big deal
CFLAGS = -O3 -fPIC -fomit-frame-pointer -funroll-loops -m64 -pedantic -std=gnu11
LDFLAGS =
INCLUDES = -I/usr/local/include
LIBS = -L/usr/local/lib -lgmp -lm
INSTALL_ROOT = /usr/local

MPZFFTHEADERS = zzcrt.h zzmem.h zzmisc.h mpzfft_moduli.h mpnfft.h mpnfft_mod.h fermat.h split.h reduce.h split_reduce.h crt.h recompose.h crt_recompose.h  fft62/arith128.h fft62/mod62.h fft62/fft62.h
MPZFFTOBJECTS = zzmisc.o moduli.o split.o reduce.o split_reduce.o crt.o recompose.o crt_recompose.o mpnfft.o fermat.o mpnfft_mod.o mpzfft.o fft62/mod62.o fft62/fft62.o zzmem.o
RFORESTHEADERS = hwmpz.h hwmpz_tune.h hwmem.h rtree.h
RFORESTOBJECTS = hwmpz.o hwmpz_tune.o hwmem.o rtree.o rforest.o
HEADERS = $(MPZFFTHEADERS) $(RFORESTHEADERS)
OBJECTS = $(MPZFFTOBJECTS) $(RFORESTOBJECTS)
PROGRAMS = test_rforest

all: librforest.a $(PROGRAMS)

clean:
	rm -f *.o
	rm -f fft62/*.o
	rm -f librforest.a $(PROGRAMS)

install: all
	cp -v rforest.h $(INSTALL_ROOT)/include
	cp -v librforest.a $(INSTALL_ROOT)/lib

##### rforest library

librforest.a: $(OBJECTS)
	ar -r librforest.a $(OBJECTS)
	ranlib librforest.a
	
##### executables

test_rforest: test_rforest.o librforest.a rforest.h
	$(CC) $(LDFLAGS) -o $@ $< librforest.a $(LIBS)

##### hwlpoly modules

hwmem.o: hwmem.c hwmem.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

hwmpz.o: hwmpz.c hwmpz.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

hwmpz_tune.o: hwmpz_tune.c hwmem.h hwmpz.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

rtree.o: rtree.c rtree.h hwmem.h hwmpz.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

rforest.o: rforest.c rforest.h hwmem.h hwmpz.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

test_rforest.o: test_rforest.c rforest.h hwmem.h hwmpz.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

##### mpzfft C modules

zzmisc.o: zzmisc.c mpzfft.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

moduli.o : moduli.c mpzfft.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

split.o : split.c mpzfft.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

reduce.o : reduce.c mpzfft.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

split_reduce.o : split_reduce.c mpzfft.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

recompose.o : recompose.c mpzfft.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

crt.o : crt.c mpzfft.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

crt_recompose.o : crt_recompose.c mpzfft.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

mpnfft.o : mpnfft.c mpzfft.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

fermat.o : fermat.c mpzfft.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

mpnfft_mod.o : mpnfft_mod.c mpzfft.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

mpzfft.o : mpzfft.c mpzfft.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

zzmem.o : zzmem.c mpzfft.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

fft62/mod62.o : fft62/mod62.c mpzfft.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

fft62/fft62.o : fft62/fft62.c mpzfft.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

##### master header files

mpzfft.h: $(MPZFFTHEADERS)
	touch mpzfft.h

rforest.h: $(HEADERS)
	touch rforest.h
